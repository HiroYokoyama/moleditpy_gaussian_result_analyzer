"""Charges, dipole, spin, dispersion, post-HF energies and NBO from a Gaussian log."""

import re
import logging

#: Debye per atomic unit of dipole moment
AU_TO_DEBYE = 2.541746473

_ATOM_CHARGE = re.compile(
    r"^\s*(\d+)\s+([A-Za-z]{1,3})\s+(-?\d+\.\d+)(?:\s+(-?\d+\.\d+))?"
)
_DIPOLE = re.compile(
    r"X=\s*(-?\d+\.\d+)\s+Y=\s*(-?\d+\.\d+)\s+Z=\s*(-?\d+\.\d+)\s+Tot=\s*(-?\d+\.\d+)"
)
_S2 = re.compile(r"<S\*\*2>=\s*(-?\d+\.\d+)")
_DISPERSION = re.compile(r"Dispersion energy=\s*(-?\d+\.\d+)\s+Hartrees")
_NBO_ROW = re.compile(
    r"^\s*(\d+)\.\s+((?:BD|CR|LP|LV|RY|3C|3Cn)\*?\s*\(\s*\d+\))\s*(.+?)\s+"
    r"(\d+\.\d{5})\s+(-?\d+\.\d{5})"
)
_E2_ROW = re.compile(
    r"^\s*(\d+)\.\s+(.+?)\s*/\s*(\d+)\.\s+(.+?)\s+"
    r"(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s*$"
)
_D_FLOAT = r"(-?\d*\.\d+(?:[DdEe][-+]?\d+)?)"


def _float_d(text):
    """Float from Fortran notation, which writes the exponent with a D."""
    return float(text.replace("D", "E").replace("d", "e"))


class _PropertyParsingMixin:
    """Parsing of scalar and per-atom properties."""

    def _last_index(self, predicate):
        """Index of the last line satisfying *predicate*, or None."""
        for i in range(len(self.lines) - 1, -1, -1):
            if predicate(self.lines[i]):
                return i
        return None

    def parse_dipole(self):
        """The last dipole moment, printed in Debye in the standard orientation."""
        self.data["dipoles"] = None
        self.data["dipole"] = None
        i = self._last_index(lambda l: "Dipole moment (field-independent basis" in l)
        if i is None or i + 1 >= len(self.lines):
            return
        match = _DIPOLE.search(self.lines[i + 1])
        if not match:
            return
        x, y, z, total = (float(match.group(k)) for k in (1, 2, 3, 4))
        self.data["dipoles"] = {
            "vector_au": (x / AU_TO_DEBYE, y / AU_TO_DEBYE, z / AU_TO_DEBYE),
            "vector_debye": (x, y, z),
            "magnitude_au": total / AU_TO_DEBYE,
            "magnitude_debye": total,
            "vector": (x, y, z),
            "magnitude": total,
        }
        self.data["dipole"] = self.data["dipoles"]

    def parse_spin_contamination(self):
        """<S**2> of the final wavefunction against the ideal S(S+1)."""
        self.data["spin_s2"] = None
        actual = None
        for line in self.lines:
            match = _S2.search(line)
            if match and "Excited State" not in line:
                actual = float(match.group(1))
        if actual is None:
            return
        s = (self.data.get("mult", 1) - 1) / 2.0
        ideal = s * (s + 1)
        self.data["spin_s2"] = {
            "actual": actual,
            "ideal": ideal,
            "contamination": actual - ideal,
        }

    def parse_dispersion(self):
        """The empirical dispersion energy, when one was requested."""
        self.data["dispersion"] = None
        for line in self.lines:
            match = _DISPERSION.search(line)
            if match:
                self.data["dispersion"] = float(match.group(1))

    def parse_energy_components(self):
        """Post-HF energies, in the order the methods build on each other."""
        found = {}
        patterns = (
            ("MP2 Correlation", re.compile(r"\bE2\s*=\s*" + _D_FLOAT), False),
            ("MP2 Total", re.compile(r"EUMP2\s*=\s*" + _D_FLOAT), False),
            (
                "CCSD Total",
                re.compile(r"Wavefunction amplitudes converged\. E\(Corr\)=\s*" + _D_FLOAT),
                False,
            ),
            ("CCSD(T) Total", re.compile(r"CCSD\(T\)=\s*" + _D_FLOAT), False),
            ("T1 Diagnostic", re.compile(r"T1 Diagnostic\s*=\s*" + _D_FLOAT), True),
        )
        for line in self.lines:
            for label, pattern, dimensionless in patterns:
                match = pattern.search(line)
                if match:
                    found[label] = (_float_d(match.group(1)), dimensionless)
        self.data["energy_components"] = [
            {"label": label, "value": found[label][0], "dimensionless": found[label][1]}
            for label, _p, _d in patterns
            if label in found
        ]

    def parse_nbo_orbitals(self):
        """The Natural Bond Orbitals summary, every spin section included."""
        orbitals = []
        for i, line in enumerate(self.lines):
            if "Natural Bond Orbitals (Summary):" not in line:
                continue
            for row in self.lines[i + 1 : i + 5000]:
                if row.strip().startswith("Total Lewis"):
                    break
                match = _NBO_ROW.match(row)
                if match:
                    atoms = match.group(3).strip()
                    orbitals.append(
                        {
                            "index": int(match.group(1)),
                            "type": " ".join(match.group(2).split()),
                            "atoms": atoms,
                            "atom_indices": [
                                int(n) - 1 for n in re.findall(r"\b\d+\b", atoms)
                            ],
                            "occupancy": float(match.group(4)),
                            "energy": float(match.group(5)),
                        }
                    )
        self.data["nbo_orbitals"] = orbitals

    def parse_nbo_perturbation(self):
        """Donor/acceptor rows of the second-order perturbation table."""
        interactions = []
        for i, line in enumerate(self.lines):
            if "Second Order Perturbation Theory Analysis" not in line:
                continue
            for row in self.lines[i + 1 : i + 20000]:
                if "Natural Bond Orbitals (Summary)" in row:
                    break
                match = _E2_ROW.match(row)
                if match:
                    interactions.append(
                        {
                            "donor": f"{match.group(1)}. {' '.join(match.group(2).split())}",
                            "acceptor": f"{match.group(3)}. {' '.join(match.group(4).split())}",
                            "e2_kcal": float(match.group(5)),
                            "e_diff": float(match.group(6)),
                            "fock": float(match.group(7)),
                        }
                    )
        self.data["nbo_perturbation"] = interactions

    def _charge_rows(self, start, spin_column=False):
        """Rows of an "atom symbol value [spin]" table beginning at *start*."""
        rows = []
        for line in self.lines[start:]:
            match = _ATOM_CHARGE.match(line)
            if not match:
                if rows:
                    break
                if line.strip() and not re.match(r"^\s*\d+(\s+\d+)*\s*$", line):
                    break
                continue
            row = {
                "atom_idx": int(match.group(1)) - 1,
                "atom_sym": match.group(2),
                "charge": float(match.group(3)),
            }
            if spin_column and match.group(4) is not None:
                row["spin"] = float(match.group(4))
            rows.append(row)
        return rows

    def _esp_scheme(self):
        """Display name of the electrostatic-potential charges this job requested."""
        route = " ".join(line for line in self.lines[:400] if line.startswith(" #"))
        route = route.upper()
        if "CHELPG" in route:
            return "CHELPG"
        if re.search(r"\bMK\b", route):
            return "Merz-Kollman"
        return "ESP"

    def parse_charges(self):
        """Every charge scheme the job printed; the last print of each wins."""
        self.data["charges"] = {}
        charges = self.data["charges"]
        esp_name = self._esp_scheme()
        npa_read = False
        for i, line in enumerate(self.lines):
            text = line.strip()
            if text in ("Mulliken charges:", "Mulliken charges and spin densities:"):
                rows = self._charge_rows(i + 1, spin_column=True)
                if rows:
                    charges["Mulliken"] = rows
            elif text == "ESP charges:":
                rows = self._charge_rows(i + 1)
                if rows:
                    charges[esp_name] = rows
            elif text == "APT charges:":
                rows = self._charge_rows(i + 1)
                if rows:
                    charges["APT"] = rows
            elif text.startswith("Hirshfeld charges, spin densities, dipoles, and CM5"):
                self._read_hirshfeld(i + 2)
            elif text.startswith("SCF Done:"):
                npa_read = False
            elif text.startswith("Summary of Natural Population Analysis") and not npa_read:
                # An open-shell job repeats this table per spin afterwards, with
                # spin populations rather than charges: only the first is charges.
                npa_read = self._read_natural_population(i)

    def _read_hirshfeld(self, start):
        """Hirshfeld and CM5 columns of the one table that prints them both."""
        hirshfeld, cm5 = [], []
        for line in self.lines[start:]:
            parts = line.split()
            if len(parts) != 8 or not parts[0].isdigit():
                break
            try:
                base = {"atom_idx": int(parts[0]) - 1, "atom_sym": parts[1]}
                hirshfeld.append(dict(base, charge=float(parts[2]), spin=float(parts[3])))
                cm5.append(dict(base, charge=float(parts[7])))
            except ValueError as exc:
                logging.warning("Hirshfeld: bad row %r: %s", line, exc)
                return
        if hirshfeld:
            self.data["charges"]["Hirshfeld"] = hirshfeld
            self.data["charges"]["CM5"] = cm5

    def _read_natural_population(self, start):
        """NPA charges from the table at *start*; True when any were read."""
        rows = []
        for line in self.lines[start + 1 : start + 20000]:
            if "Total *" in line:
                break
            parts = line.split()
            if len(parts) >= 7 and parts[1].isdigit():
                try:
                    rows.append(
                        {
                            "atom_idx": int(parts[1]) - 1,
                            "atom_sym": parts[0],
                            "charge": float(parts[2]),
                        }
                    )
                except ValueError as exc:
                    logging.warning("NPA: bad row %r: %s", line, exc)
                    return False
        if rows:
            self.data["charges"]["NBO"] = rows
        return bool(rows)
