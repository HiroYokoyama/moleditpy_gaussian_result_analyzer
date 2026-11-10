"""Geometry, job status, forces and optimization/scan steps from a Gaussian log."""

from __future__ import annotations

import re
import math

try:
    import numpy as np
except ImportError:  # the pytest-only CI job has no numpy
    np = None

#: Atomic symbols indexed by atomic number; 0 is a ghost/dummy centre.
ELEMENTS = (
    "X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni "
    "Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe "
    "Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au "
    "Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf "
    "Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og"
).split()

_GEOMETRY_HEADER = re.compile(r"^\s*(Input|Standard|Z-Matrix) orientation:\s*$")
_STEP = re.compile(
    r"^\s*Step number\s+(\d+)\s+out of a maximum of\s+\d+"
    r"(?:\s+on scan point\s+(\d+)\s+out of\s+(\d+))?"
)
_SCF_DONE = re.compile(r"^\s*SCF Done:\s+E\(\S+\)\s*=\s*(-?\d+\.\d+)")
_CONVERGED_ROW = re.compile(
    r"^\s*((?:Maximum|RMS)\s+(?:Force|Displacement))\s+"
    r"(-?\d+\.\d+|\*+)\s+(\d+\.\d+)\s+(YES|NO)\s*$"
)
_JOB_STEP = "Proceeding to internal job step"


def symbol_for(atomic_number):
    """Element symbol for an atomic number; ghosts and dummies are "X"."""
    if 0 < atomic_number < len(ELEMENTS):
        return ELEMENTS[atomic_number]
    return "X"


class _StructureParsingMixin:
    """Parsing of geometry, job status, forces and optimization trajectories."""

    def _read_geometry(self, header_idx):
        """(atoms, coords) of the orientation table whose title is at header_idx."""
        atoms, coords = [], []
        for line in self.lines[header_idx + 5 :]:
            parts = line.split()
            if len(parts) != 6:
                break
            try:
                atomic_number = int(parts[1])
                xyz = (float(parts[3]), float(parts[4]), float(parts[5]))
            except ValueError:
                break
            if not all(math.isfinite(v) for v in xyz):
                break
            atoms.append(symbol_for(atomic_number))
            coords.append(list(xyz))
        return atoms, coords

    def parse_basic(self):
        """Version, charge/multiplicity, final energy and geometry, job status."""
        last_standard = last_any = None
        for i, line in enumerate(self.lines):
            match = _GEOMETRY_HEADER.match(line)
            if match:
                last_any = i
                if match.group(1) == "Standard":
                    last_standard = i
                continue
            version = re.match(r"^ (Gaussian \d+):\s+(\S+)", line)
            if version:
                self.data["version"] = f"{version.group(1)} {version.group(2)}"
                continue
            charge = re.match(
                r"^\s*Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)", line
            )
            if charge:
                self.data["charge"] = int(charge.group(1))
                self.data["mult"] = int(charge.group(2))
                continue
            scf = _SCF_DONE.match(line)
            if scf:
                self.data["scf_energy"] = float(scf.group(1))
                continue
            if "Stationary point found" in line or "Optimization completed." in line:
                self.data["converged"] = True
            if "on scan point" in line:
                self.data["is_scan"] = True

        # MO coefficients, normal modes and the dipole are all in the standard
        # orientation; it is only absent under NoSymm, where input is used.
        header = last_standard if last_standard is not None else last_any
        if header is not None:
            atoms, coords = self._read_geometry(header)
            if atoms:
                self.data["atoms"], self.data["coords"] = atoms, coords

        self.data["termination_status"] = "Running"
        for line in reversed(self.lines[-400:]):
            if "Normal termination of Gaussian" in line:
                self.data["termination_status"] = "Terminated normally"
                break
            if "Error termination" in line:
                self.data["termination_status"] = "ERROR"
                break

    @staticmethod
    def _rotation_between(source, target):
        """Rotation taking source coordinates onto target (Kabsch), or None."""
        if np is None or len(source) != len(target) or len(source) < 2:
            return None
        a = np.array(source, dtype=float)
        b = np.array(target, dtype=float)
        a = a - a.mean(axis=0)
        b = b - b.mean(axis=0)
        u, _s, vt = np.linalg.svd(a.T @ b)
        d = np.sign(np.linalg.det(vt.T @ u.T))
        rotation = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
        if np.abs((rotation @ a.T).T - b).max() > 0.05:
            return None  # not the same structure: do not invent a frame
        return rotation

    def _input_section_coords(self):
        """Cartesian coordinates as typed in the input, or None for a Z-matrix."""
        start = next(
            (i for i, l in enumerate(self.lines) if "Symbolic Z-matrix:" in l), None
        )
        if start is None:
            return None
        coords = []
        for line in self.lines[start + 2 :]:
            parts = line.split()
            if not parts:
                break
            if len(parts) == 5 and parts[1] in ("0", "-1"):
                parts.pop(1)  # the freeze flag
            try:
                coords.append([float(v) for v in parts[1:4]])
            except ValueError:
                return None
            if len(parts) != 4:
                return None
        return coords or None

    def parse_gradients(self):
        """Read every force block as gradients (schema: gradient = -force).

        Gaussian prints forces in the input orientation, while data["coords"]
        is in the standard orientation, so the final gradient is rotated onto
        it. Otherwise the arrows would point the wrong way on the molecule.
        """
        self.data["gradients"] = []
        self.data["all_gradients"] = []
        last_input = None
        for i, line in enumerate(self.lines):
            match = _GEOMETRY_HEADER.match(line)
            if match and match.group(1) in ("Input", "Z-Matrix"):
                last_input = i
            if "Forces (Hartrees/Bohr)" not in line:
                continue
            block = self._read_force_block(i)
            if block:
                self.data["all_gradients"].append(
                    {"line": i, "grads": block, "geometry_line": last_input}
                )
        if not self.data["all_gradients"]:
            return
        final = self.data["all_gradients"][-1]
        grads = [dict(g, vector=list(g["vector"])) for g in final["grads"]]
        input_coords = None
        if final["geometry_line"] is not None:
            _atoms, input_coords = self._read_geometry(final["geometry_line"])
        elif len(self.data["all_gradients"]) == 1:
            # Above 50 atoms the input orientation is not printed.
            input_coords = self._input_section_coords()
        if input_coords and self.data.get("coords"):
            rotation = self._rotation_between(input_coords, self.data["coords"])
            if rotation is not None:
                for g in grads:
                    g["vector"] = [float(v) for v in rotation @ np.array(g["vector"])]
        self.data["gradients"] = grads

    def _read_force_block(self, header_idx):
        """The gradient rows of the force table titled at header_idx."""
        block = []
        for line in self.lines[header_idx + 3 :]:
            parts = line.split()
            if len(parts) != 5:
                break
            try:
                idx, atomic_number = int(parts[0]), int(parts[1])
                force = [float(parts[2]), float(parts[3]), float(parts[4])]
            except ValueError:
                break
            block.append(
                {
                    "atom_idx": idx - 1,
                    "atom_sym": symbol_for(atomic_number),
                    "vector": [-v for v in force],
                }
            )
        return block

    def _job_ranges(self):
        """(start, end) line ranges of each Link1 job step."""
        starts = [0] + [i for i, l in enumerate(self.lines) if _JOB_STEP in l]
        return list(zip(starts, starts[1:] + [len(self.lines)]))

    def parse_trajectory(self):
        """One entry per optimization point in data["scan_steps"].

        A job step without a "Step number" line is not an optimization: the
        Freq step of Opt Freq also prints forces and a convergence table, and
        must not appear as one more cycle. Coordinates stay in the input
        orientation, the frame the forces are printed in; the standard
        orientation may be rotated from one point to the next.
        """
        self.data["scan_steps"] = []
        for start, end in self._job_ranges():
            if any(_STEP.match(self.lines[i]) for i in range(start, end)):
                self._read_optimization(start, end)

    def _read_optimization(self, start, end):
        """Append the optimization points found in lines[start:end]."""
        steps = self.data["scan_steps"]
        grads_at = {g["line"]: g["grads"] for g in self.data["all_gradients"]}
        current = {"energy": None, "geometry": None, "grads": [], "conv": {}}
        counts = {}
        scan = bool(self.data.get("is_scan"))

        def flush(point):
            if current["geometry"] is None or current["energy"] is None:
                return
            atoms, coords = current["geometry"]
            counts[point] = counts.get(point, 0) + 1
            steps.append(
                {
                    "type": "scan_step" if scan else "opt_cycle",
                    "scan_step_id": point,
                    "step": counts[point],
                    "energy": current["energy"],
                    "atoms": atoms,
                    "coords": coords,
                    "convergence": dict(current["conv"]),
                    "gradients": current["grads"],
                }
            )
            current.update(energy=None, geometry=None, grads=[], conv={})

        point = 1
        for i in range(start, end):
            line = self.lines[i]
            header = _GEOMETRY_HEADER.match(line)
            if header:
                if header.group(1) != "Standard" or current["geometry"] is None:
                    current["geometry"] = self._read_geometry(i)
                continue
            scf = _SCF_DONE.match(line)
            if scf:
                current["energy"] = float(scf.group(1))
                continue
            if i in grads_at:
                current["grads"] = grads_at[i]
                continue
            row = _CONVERGED_ROW.match(line)
            if row:
                current["conv"][" ".join(row.group(1).split())] = {
                    "value": row.group(2),
                    "tolerance": row.group(3),
                    "converged": row.group(4),
                }
                continue
            step = _STEP.match(line)
            if step:
                # This line closes the geometry above it; its scan point number
                # names that same point, so it is read before flushing.
                if step.group(2):
                    point = int(step.group(2))
                flush(point)
        flush(point)

    def parse_scan_results_table(self):
        """Attach the scanned coordinate from the end-of-job summary table."""
        scan_names = []
        for line in self.lines:
            match = re.match(r"^\s*!\s+([RADL]\d+)\s+\S+.*\bScan\b", line)
            if match and match.group(1) not in scan_names:
                scan_names.append(match.group(1))
        if not scan_names:
            return
        values = []
        for i, line in enumerate(self.lines):
            if "Summary of Optimized Potential Surface Scan" not in line:
                continue
            for row_line in self.lines[i + 1 :]:
                row = row_line.split()
                if not row or row_line.strip().startswith("---"):
                    break
                if row[0] == scan_names[0]:
                    values.extend(float(v) for v in row[1:])
            break
        for step in self.data["scan_steps"]:
            idx = step.get("scan_step_id", 0) - 1
            if 0 <= idx < len(values):
                step["scan_coord"] = values[idx]
