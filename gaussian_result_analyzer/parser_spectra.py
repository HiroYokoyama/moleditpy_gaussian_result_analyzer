"""Vibrational modes, thermochemistry, TD-DFT and NMR from a Gaussian log."""

import re
import logging

#: Below this magnitude an imaginary frequency is taken as numerical noise.
IMAGINARY_FREQ_THRESHOLD = 10.0

KCAL_PER_HARTREE = 627.5094740631
EV_TO_CM = 8065.543937

_EXCITED = re.compile(
    r"^\s*Excited State\s+(\d+):\s+(\S+)\s+(-?\d+\.\d+)\s+eV\s+(-?\d+\.\d+)\s+nm\s+"
    r"f=\s*(-?\d+\.\d+)(?:\s+<S\*\*2>=\s*(-?\d+\.\d+))?"
)
_TRANSITION = re.compile(r"^\s*(\d+[AB]?)\s*(->|<-)\s*(\d+[AB]?)\s+(-?\d+\.\d+)\s*$")
_SHIELDING = re.compile(
    r"^\s*(\d+)\s+([A-Za-z]{1,3})\s+Isotropic\s*=\s*(-?\d+\.\d+)"
    r"\s+Anisotropy\s*=\s*(-?\d+\.\d+)"
)
_D = r"[-+]?\d*\.\d+[DdEe][-+]?\d+"


def _float_d(text):
    """Float from Fortran notation, which writes the exponent with a D."""
    return float(text.replace("D", "E").replace("d", "e"))


class _SpectraParsingMixin:
    """Parsing of spectra and thermochemistry."""

    def parse_frequencies(self):
        """Normal modes from the last Harmonic frequencies section.

        Only the 3-per-row table is read. Freq=HPModes prints a higher-precision
        table first and this one after it, so the two never need merging.
        """
        self.data["frequencies"] = []
        start = None
        for i, line in enumerate(self.lines):
            if "Harmonic frequencies (cm**-1)" in line:
                start = i
        if start is None:
            return
        n_atoms = len(self.data.get("atoms", []))
        modes = []
        current = []
        i = start
        while i < len(self.lines):
            stripped = self.lines[i].strip()
            if stripped.startswith("Frequencies --"):
                values = [float(v) for v in stripped.split("--", 1)[1].split()]
                current = [
                    {"freq": v, "ir": 0.0, "raman": 0.0, "vector": []} for v in values
                ]
                modes.extend(current)
            elif stripped.startswith("Red. masses --"):
                self._column(current, "reduced_mass", stripped)
            elif stripped.startswith("Frc consts  --"):
                self._column(current, "force_constant", stripped)
            elif stripped.startswith("IR Inten    --"):
                self._column(current, "ir", stripped)
            elif stripped.startswith("Raman Activ --"):
                self._column(current, "raman", stripped)
            elif stripped.startswith("Atom  AN") and current:
                i = self._read_mode_vectors(i + 1, current, n_atoms)
                continue
            elif stripped.startswith("- Thermochemistry -"):
                break
            i += 1
        self.data["frequencies"] = modes

    @staticmethod
    def _column(modes, key, stripped):
        """Assign one numeric column of a Frequencies block to each of its modes."""
        values = stripped.split("--", 1)[1].split()
        for mode, text in zip(modes, values):
            try:
                mode[key] = float(text)
            except ValueError:
                logging.debug("Unparseable %s value %r", key, text)

    def _read_mode_vectors(self, start, modes, n_atoms):
        """Read one block's displacement rows; returns the next line index."""
        i = start
        while i < len(self.lines) and (n_atoms == 0 or i - start < n_atoms):
            parts = self.lines[i].split()
            if len(parts) != 2 + 3 * len(modes) or not parts[0].isdigit():
                break
            try:
                values = [float(v) for v in parts[2:]]
            except ValueError:
                break
            for k, mode in enumerate(modes):
                mode["vector"].append(tuple(values[k * 3 : k * 3 + 3]))
            i += 1
        return i

    def parse_thermal(self):
        """Thermochemistry, in Hartree (Gaussian prints kcal/mol and cal/mol-K)."""
        self.data["thermal"] = {}
        start = None
        for i, line in enumerate(self.lines):
            if line.strip() == "- Thermochemistry -":
                start = i
        if start is None:
            return
        data = self.data["thermal"]
        values = {}
        contributions = {}
        for line in self.lines[start:]:
            match = re.match(r"^\s*Temperature\s+(\d+\.\d+)\s+Kelvin", line)
            if match:
                data["temperature"] = float(match.group(1))
                continue
            match = re.match(
                r"^\s*(Zero-point correction|Thermal correction to [\w ]+?|"
                r"Sum of electronic and [\w -]+?)\s*=\s*(-?\d+\.\d+)",
                line,
            )
            if match:
                values[match.group(1)] = float(match.group(2))
                continue
            match = re.match(
                r"^\s*(Total|Electronic|Translational|Rotational|Vibrational)\s+"
                r"(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$",
                line,
            )
            if match and match.group(1) not in contributions:
                contributions[match.group(1)] = (
                    float(match.group(2)),
                    float(match.group(4)),
                )
            if "Total Bot" in line:
                break
        zpe = values.get("Zero-point correction")
        zpe_sum = values.get("Sum of electronic and zero-point Energies")
        if zpe is None or zpe_sum is None:
            return
        data["zpe"] = zpe
        data["electronic_energy"] = zpe_sum - zpe
        data["thermal_energy"] = values.get("Sum of electronic and thermal Energies")
        data["enthalpy"] = values.get("Sum of electronic and thermal Enthalpies")
        data["gibbs"] = values.get("Sum of electronic and thermal Free Energies")
        data["enthalpy_corr"] = values.get("Thermal correction to Enthalpy")
        data["gibbs_corr"] = values.get("Thermal correction to Gibbs Free Energy")
        energy_corr = values.get("Thermal correction to Energy")
        temperature = data.get("temperature", 298.15)

        def ts(key):
            """T*S in Hartree for one row of the contribution table."""
            if key not in contributions:
                return None
            return temperature * contributions[key][1] / 1000.0 / KCAL_PER_HARTREE

        def thermal(key):
            """E(thermal) in Hartree for one row of the contribution table."""
            if key not in contributions:
                return None
            return contributions[key][0] / KCAL_PER_HARTREE

        data["entropy"] = ts("Total")
        data["s_el"] = ts("Electronic")
        data["s_vib"] = ts("Vibrational")
        data["s_rot"] = ts("Rotational")
        data["s_trans"] = ts("Translational")
        vib = thermal("Vibrational")
        # Gaussian's vibrational E(thermal) includes the zero-point energy; the
        # table this fills reports the two separately.
        data["corr_vib"] = vib - zpe if vib is not None else None
        data["corr_rot"] = thermal("Rotational")
        data["corr_trans"] = thermal("Translational")
        data["corr_zpe"] = zpe
        data["corr_total"] = energy_corr
        if energy_corr is not None:
            data["corr_thermal_total"] = energy_corr - zpe
        data["imaginary_freq_count"] = sum(
            1
            for f in self.data.get("frequencies", [])
            if f.get("freq", 0) < -IMAGINARY_FREQ_THRESHOLD
        )

    def parse_tddft(self):
        """Excited states of the last TD/CIS section, length and velocity gauge."""
        self.data["tddft"] = []
        states = {}
        order = []
        start = None
        for i, line in enumerate(self.lines):
            if "Excitation energies and oscillator strengths:" in line:
                start = i
        if start is None:
            return
        current = None
        for line in self.lines[start:]:
            match = _EXCITED.match(line)
            if match:
                ev = float(match.group(3))
                state = {
                    "state": int(match.group(1)),
                    "label": match.group(2),
                    "energy_ev": ev,
                    "energy_nm": float(match.group(4)),
                    "energy_cm": ev * EV_TO_CM,
                    "osc": float(match.group(5)),
                    "osc_len": float(match.group(5)),
                    "osc_vel": 0.0,
                    "rotatory_strength": 0.0,
                    "rot_len": 0.0,
                    "rot_vel": 0.0,
                    "transitions": [],
                }
                if match.group(6) is not None:
                    state["s2"] = float(match.group(6))
                states[state["state"]] = state
                order.append(state["state"])
                current = state
                continue
            if current is not None:
                match = _TRANSITION.match(line)
                if match:
                    current["transitions"].append(
                        f"{match.group(1)} {match.group(2)} {match.group(3)} "
                        f"(coeff: {match.group(4)})"
                    )
        self._read_gauge_tables(states, start)
        self.data["tddft"] = [states[k] for k in order]

    def _read_gauge_tables(self, states, before):
        """Velocity-gauge and rotatory strengths from the tables printed above
        the excited state list. Each table announces its columns in a "state"
        header, which is what decides how its rows are read."""
        begin = max(0, before - 12 * (len(states) + 4))
        section = None
        for line in self.lines[begin:before]:
            parts = line.split()
            if parts and parts[0] == "state":
                if "Osc.(frdel)" in line:
                    section = "vel"
                elif "R(length)" in line:
                    section = "rot_len"
                elif "R(velocity)" in line:
                    section = "rot_vel"
                else:
                    section = None
                continue
            if section is None or len(parts) < 5 or not parts[0].isdigit():
                continue
            state = states.get(int(parts[0]))
            if state is None:
                continue
            try:
                if section == "vel":
                    state["osc_vel"] = float(parts[-1])
                elif section == "rot_vel":
                    state["rot_vel"] = float(parts[4])
                else:
                    state["rot_len"] = float(parts[4])
                    state["rotatory_strength"] = state["rot_len"]
            except ValueError:
                logging.debug("TD: unparseable table row %r", line)

    def parse_nmr(self):
        """Isotropic shieldings and, with NMR=SpinSpin, total J couplings in Hz."""
        self.data["nmr_shielding"] = []
        self.data["nmr_couplings"] = []
        start = None
        for i, line in enumerate(self.lines):
            if "Magnetic shielding tensor (ppm)" in line:
                start = i
        if start is None:
            return
        for line in self.lines[start:]:
            match = _SHIELDING.match(line)
            if match:
                self.data["nmr_shielding"].append(
                    {
                        "atom_idx": int(match.group(1)) - 1,
                        "atom_sym": match.group(2),
                        "shielding": float(match.group(3)),
                        "anisotropy": float(match.group(4)),
                    }
                )
        for i, line in enumerate(self.lines):
            if line.strip().startswith("Total nuclear spin-spin coupling J (Hz)"):
                self._read_lower_triangle(i + 1)

    def _read_lower_triangle(self, start):
        """Read a lower-triangular J matrix, printed in blocks of five columns."""
        couplings = {}
        columns = []
        for line in self.lines[start:]:
            parts = line.split()
            if not parts:
                break
            if all(p.isdigit() for p in parts):
                columns = [int(p) - 1 for p in parts]
                continue
            values = re.findall(_D, line)
            if not parts[0].isdigit() or not values or not columns:
                break
            row = int(parts[0]) - 1
            for column, text in zip(columns, values):
                if column < row:
                    couplings[(column, row)] = _float_d(text)
        self.data["nmr_couplings"] = [
            {"atom_idx1": a, "atom_idx2": b, "coupling": j}
            for (a, b), j in sorted(couplings.items())
        ]
