"""MO coefficients, orbital energies, basis set and SCF traces from a Gaussian log."""

import re
import logging

HARTREE_TO_EV = 27.211386245988

_SHELL_L = {"S": 0, "P": 1, "D": 2, "F": 3, "G": 4, "H": 5, "I": 6}
_MO_NUMBER = re.compile(r"-?\d*\.\d{5}")
_NUMBER_D = re.compile(r"[-+]?\d*\.?\d+(?:[DdEe][-+]?\d+)?")
#: "  16 2   H  1S   ..." or "   2        2S   ...": the orbital label ends where
#: the fixed-width F10.5 columns begin, and those may touch each other.
_MO_ROW = re.compile(
    r"^\s*(\d+)\s+(?:(\d+)\s+([A-Za-z]{1,3})\s+)?"
    r"(\d+[A-Z](?:[XYZ]+|\s?[+-]?\d+)?)\s+(-?\d*\.\d{5}.*)$"
)
_ELECTRONS = re.compile(r"(\d+)\s+alpha electrons\s+(\d+)\s+beta electrons")
_CYCLE = re.compile(r"^\s*Cycle\s+(\d+)\s+Pass")
_CYCLE_ENERGY = re.compile(r"^\s*E=\s*(-?\d+\.\d+)")


def _float_d(text):
    """Float from Fortran notation, which writes the exponent with a D."""
    return float(text.replace("D", "E").replace("d", "e"))


def _is_number(text):
    """Whether *text* parses as a float (an overflowed field is all asterisks)."""
    try:
        float(text)
    except ValueError:
        return False
    return True


class _ElectronicParsingMixin:
    """Parsing of orbitals, the basis set and SCF convergence."""

    def _electron_counts(self):
        """(n_alpha, n_beta) from the last count printed, or (None, None)."""
        for line in reversed(self.lines):
            match = _ELECTRONS.search(line)
            if match:
                return int(match.group(1)), int(match.group(2))
        return None, None

    def _occupation(self, index, spin):
        """Electrons in 0-based orbital *index*.

        The printed O/V tag is not used: it does not distinguish a singly from
        a doubly occupied orbital, which matters for ROHF.
        """
        n_alpha, n_beta = self._electron_counts()
        if n_alpha is None:
            return 0.0
        if spin == "alpha":
            return 1.0 if index < n_alpha else 0.0
        if spin == "beta":
            return 1.0 if index < n_beta else 0.0
        return float(index < n_alpha) + float(index < n_beta)

    def _last_mo_blocks(self):
        """{spin: title line index}, keeping only the last print of each spin."""
        found = {}
        for i, line in enumerate(self.lines):
            text = line.strip()
            if text == "Molecular Orbital Coefficients:":
                found = {"restricted": i}
            elif text == "Alpha Molecular Orbital Coefficients:":
                found = {"alpha": i}
            elif text == "Beta Molecular Orbital Coefficients:":
                found["beta"] = i
        return found

    def parse_mo_coeffs(self) -> None:
        """Extract per-orbital coefficients, energy, occupancy and spin.

        Keys are "<0-based index>_<spin>". Coefficients are in the standard
        orientation and in the order Gaussian prints its basis functions.
        """
        self.data["mo_coeffs"] = {}
        self.data["cartesian_d"] = False
        self.data["cartesian_f"] = False
        for spin, title in self._last_mo_blocks().items():
            self._read_mo_block(title + 1, spin)

    def _read_mo_block(self, start, spin):
        """Read the column blocks of one MO table, starting at line *start*."""
        columns = []
        atom_idx = -1
        sym = ""
        for line in self.lines[start:]:
            parts = line.split()
            if not parts:
                break
            if all(p.isdigit() for p in parts):
                columns = [int(p) - 1 for p in parts]
                continue
            stripped = line.strip()
            if stripped.startswith("Eigenvalues --"):
                # Overflowing energies print as asterisks; those stay at 0.
                values = line.split("--", 1)[1]
                for k, column in enumerate(columns):
                    text = values[k * 10 : k * 10 + 10].strip()
                    self.data["mo_coeffs"][f"{column}_{spin}"] = {
                        "coeffs": [],
                        "spin": spin,
                        "energy": float(text) if _is_number(text) else 0.0,
                        "occ": self._occupation(column, spin),
                        "id": column,
                    }
                continue
            row = _MO_ROW.match(line)
            if row and columns:
                if row.group(2):
                    atom_idx = int(row.group(2)) - 1
                    sym = row.group(3)
                orb = row.group(4)
                self._note_cartesian(orb)
                for k, text in enumerate(_MO_NUMBER.findall(row.group(5))):
                    if k >= len(columns):
                        break
                    self.data["mo_coeffs"][f"{columns[k]}_{spin}"]["coeffs"].append(
                        {
                            "atom_idx": atom_idx,
                            "sym": sym,
                            "orb": orb,
                            "coeff": float(text),
                        }
                    )
                continue
            if re.match(r"^(\(\S+\)(--[OV]?)?\s*|[OV]\s*)+$", stripped):
                continue  # symmetry/occupancy; long labels like (DLTG)-- lose the O/V
            break  # anything else ends the table

    def _note_cartesian(self, orb):
        """Record whether d/f shells are Cartesian, from the function labels."""
        label = re.sub(r"^\d+", "", orb)
        if label in ("XX", "YY", "ZZ", "XY", "XZ", "YZ"):
            self.data["cartesian_d"] = True
        elif len(label) == 3 and set(label) <= set("XYZ"):
            self.data["cartesian_f"] = True

    def parse_orbital_energies(self):
        """Orbital energies, from the MO table or else the eigenvalue summary."""
        self.data["orbital_energies"] = []
        self.data["mos"] = []
        entries = []
        if self.data.get("mo_coeffs"):
            for mo in self.data["mo_coeffs"].values():
                entries.append((mo["spin"], mo["id"], mo["energy"], mo["occ"]))
        else:
            entries = self._energies_from_summary()
        for spin, index, energy, occ in sorted(entries, key=lambda e: (e[0], e[1])):
            orb = {
                "index": index,
                "id": index,
                "occupation": occ,
                "occ": occ,
                "energy_eh": energy,
                "energy_ev": energy * HARTREE_TO_EV,
                "energy": energy,
                "spin": spin,
                "type": "occupied" if occ > 0.1 else "virtual",
            }
            self.data["orbital_energies"].append(orb)
            self.data["mos"].append(orb)

    def _energies_from_summary(self):
        """(spin, index, energy, occ) from the last "occ./virt. eigenvalues" lines."""
        last_block = None
        for i, line in enumerate(self.lines):
            if "occ. eigenvalues --" in line and (
                i == 0 or "eigenvalues --" not in self.lines[i - 1]
            ):
                last_block = i
        if last_block is None:
            return []
        values = {"alpha": [], "beta": []}
        for line in self.lines[last_block:]:
            if "eigenvalues --" not in line:
                break
            spin = "beta" if line.strip().startswith("Beta") else "alpha"
            tail = line.split("--", 1)[1]
            for k in range(0, len(tail), 10):
                text = tail[k : k + 10].strip()
                if _is_number(text):
                    values[spin].append(float(text))
        out = []
        names = ["alpha", "beta"] if values["beta"] else ["restricted"]
        for name in names:
            source = values["alpha" if name == "restricted" else name]
            out.extend(
                (name, i, e, self._occupation(i, name)) for i, e in enumerate(source)
            )
        return out

    def parse_basis_set(self):
        """Contracted shells from the GFInput block.

        An SP shell becomes an S and a P shell sharing exponents, which is the
        order Gaussian numbers those basis functions in.
        """
        self.data["basis_set_shells"] = []
        start = None
        for i, line in enumerate(self.lines):
            if "AO basis set in the form of general basis input" in line:
                start = i
        coords = self.data.get("coords", [])
        if start is None or not coords:
            return
        shells = []
        atom_idx = None
        i = start + 1
        while i < len(self.lines):
            parts = self.lines[i].split()
            if len(parts) == 2 and parts[1] == "0" and parts[0].isdigit():
                atom_idx = int(parts[0]) - 1
                i += 1
                continue
            if parts == ["****"]:
                i += 1
                continue
            if len(parts) < 3 or not re.fullmatch(r"[SPDFGHI]+", parts[0]):
                break
            try:
                n_prim = int(parts[1])
                scale = _float_d(parts[2])
                rows = [
                    [_float_d(v) for v in _NUMBER_D.findall(self.lines[i + 1 + k])]
                    for k in range(n_prim)
                ]
                exps = [row[0] * scale * scale for row in rows]
                columns = [
                    [row[column + 1] for row in rows]
                    for column in range(len(parts[0]))
                ]
            except (IndexError, ValueError) as exc:
                logging.warning("Basis set: bad shell at line %d: %s", i, exc)
                return
            if atom_idx is None or atom_idx >= len(coords):
                logging.warning("Basis set: shell for an unknown atom at line %d", i)
                return
            for letter, coeffs in zip(parts[0], columns):
                cartesian = (letter == "D" and self.data.get("cartesian_d")) or (
                    letter == "F" and self.data.get("cartesian_f")
                )
                shells.append(
                    {
                        "atom_idx": atom_idx,
                        "origin": coords[atom_idx],
                        "l": _SHELL_L[letter],
                        "exps": exps,
                        "coeffs": coeffs,
                        "cartesian": bool(cartesian),
                    }
                )
            i += 1 + n_prim
        self.data["basis_set_shells"] = shells

    def parse_scf_trace(self):
        """Energy per SCF cycle, one trace for each SCF in the job."""
        self.data["scf_traces"] = []
        trace = []
        cycle = None

        def finish():
            if trace:
                n = len(self.data["scf_traces"]) + 1
                self.data["scf_traces"].append(
                    {"step": f"SCF {n}", "iterations": list(trace)}
                )
                trace.clear()

        for line in self.lines:
            match = _CYCLE.match(line)
            if match:
                new_cycle = int(match.group(1))
                if cycle is not None and new_cycle <= cycle:
                    finish()
                cycle = new_cycle
                continue
            if cycle is not None:
                energy = _CYCLE_ENERGY.match(line)
                if energy:
                    trace.append({"iter": cycle, "energy": float(energy.group(1))})
                    continue
            if line.startswith(" SCF Done:"):
                finish()
                cycle = None
        finish()
