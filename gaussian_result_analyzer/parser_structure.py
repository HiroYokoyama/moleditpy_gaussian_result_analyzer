"""Geometry, job status, forces and optimization/scan steps from a Gaussian log."""

from __future__ import annotations

import logging
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
#: The converged energy: SCF Done, or an MCSCF (CASSCF) iteration, which prints none.
_SCF_DONE = re.compile(
    r"^\s*(?:SCF Done:\s+E\(\S+\)\s*=|ITN=\s*\d+\s+MaxIt=\s*\d+\s+E=)"
    r"\s*(-?\d+\.\d+)"
)
_CONVERGED_ROW = re.compile(
    r"^\s*((?:Maximum|RMS)\s+(?:Force|Displacement))\s+"
    r"(-?\d+\.\d+|\*+)\s+(\d+\.\d+)\s+(YES|NO)\s*$"
)
_JOB_STEP = "Proceeding to internal job step"
_VERSION = re.compile(r"^ (Gaussian \d+):\s+(\S+)")
_OPT_KEYWORD = re.compile(r"\b[fp]?opt\b", re.I)


def _ends_job(line):
    """True on the line that closes a job step; the next one may follow."""
    return _JOB_STEP in line or "Normal termination of Gaussian" in line


def _route_text(lines):
    """The first route among *lines* as one string, or None without one.

    Gaussian echoes a long route in fixed-width rows that may break mid-word.
    """
    for i, line in enumerate(lines):
        if line.startswith(" #"):
            route = []
            for row in lines[i:]:
                if row.strip().startswith("---"):
                    break
                route.append(row[1:].rstrip("\r\n"))
            return "".join(route)
    return None


_IRC_BANNER = "IRC-IRC-IRC"
_IRC_POINT = re.compile(r"^\s*Point Number:\s+(\d+)\s+Path Number:\s+(\d+)")
_IRC_DIRECTION = re.compile(r"^\s*Point Number\s+\d+\s+in\s+(FORWARD|REVERSE)\s+path")
_IRC_NET = re.compile(r"NET REACTION COORDINATE UP TO THIS POINT\s*=\s*(-?\d+\.\d+)")
_IRC_CHK_ENERGY = re.compile(r"^\s*Energy From Chk\s*=\s*(-?\d+\.\d+)")
_RIGID_SCAN = "Scan the potential surface."


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
        new_job = False
        charge_above = False
        for i, line in enumerate(self.lines):
            if _ends_job(line):
                new_job = True
            match = _GEOMETRY_HEADER.match(line)
            if match:
                if new_job:
                    # The final structure is that of the last job step printing one.
                    last_standard, new_job = None, False
                last_any = i
                if match.group(1) == "Standard":
                    last_standard = i
                continue
            charge = re.match(
                r"^\s*Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)", line
            )
            if charge:
                # Counterpoise fragments and ONIOM layers follow the whole system's.
                if not charge_above:
                    self.data["charge"] = int(charge.group(1))
                    self.data["mult"] = int(charge.group(2))
                charge_above = True
                continue
            charge_above = False
            version = _VERSION.match(line)
            if version:
                self.data["version"] = f"{version.group(1)} {version.group(2)}"
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
            if _JOB_STEP in line or _VERSION.match(line):
                break  # a later job step began and has not ended

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

    @staticmethod
    def _best_fit(coords, reference):
        """coords rigidly moved onto reference by least squares."""
        a = np.array(coords, dtype=float)
        b = np.array(reference, dtype=float)
        u, _s, vt = np.linalg.svd((a - a.mean(axis=0)).T @ (b - b.mean(axis=0)))
        d = np.sign(np.linalg.det(vt.T @ u.T))
        rotation = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
        return ((rotation @ (a - a.mean(axis=0)).T).T + b.mean(axis=0)).tolist()

    @staticmethod
    def _torque_free(coords, vectors):
        """vectors turned slightly so that their net torque about coords vanishes.

        The energy does not change when the molecule rotates, so its gradient
        has no net torque in the frame it was computed in. None when no small
        rotation achieves that.
        """
        x = np.array(coords, dtype=float)
        x = x - x.mean(axis=0)
        g = np.array(vectors, dtype=float)
        total = np.zeros(3)
        for _ in range(20):
            torque = np.cross(x, g).sum(axis=0)
            if np.abs(torque).max() < 1e-9:
                return g.tolist()
            jacobian = np.einsum("ij,ij->", x, g) * np.eye(3) - g.T @ x
            try:
                w = np.linalg.solve(jacobian, -torque)
            except np.linalg.LinAlgError:
                return None
            total += w
            if np.linalg.norm(total) > 0.3:
                return None
            angle = np.linalg.norm(w)
            if angle:
                k = w / angle
                g = (
                    g * math.cos(angle)
                    + np.cross(k, g) * math.sin(angle)
                    + np.outer(g @ k, k) * (1 - math.cos(angle))
                )
        return None

    def _synthetic_input_frames(self):
        """{standard table line: coords near the input orientation} above 50 atoms.

        Gaussian then prints only the standard orientation, while the forces
        stay in the input one. Each table is fitted onto the input section;
        the forces read beside it are then turned the last small way.
        """
        frames = {}
        if np is None or any(
            (m := _GEOMETRY_HEADER.match(l)) and m.group(1) != "Standard"
            for l in self.lines
        ):
            return frames
        reference = self._input_section_coords()
        if not reference:
            return frames
        for i, line in enumerate(self.lines):
            if _GEOMETRY_HEADER.match(line):
                _atoms, coords = self._read_geometry(i)
                if len(coords) == len(reference):
                    frames[i] = self._best_fit(coords, reference)
        return frames

    def parse_gradients(self):
        """Read every force block as gradients (schema: gradient = -force).

        Gaussian prints forces in the input orientation, while data["coords"]
        is in the standard orientation, so the final gradient is rotated onto
        it. Otherwise the arrows would point the wrong way on the molecule.
        """
        self.data["gradients"] = []
        self.data["all_gradients"] = []
        self._input_frames = self._synthetic_input_frames()
        last_input = None
        for i, line in enumerate(self.lines):
            match = _GEOMETRY_HEADER.match(line)
            if match and (match.group(1) != "Standard" or i in self._input_frames):
                last_input = i
            if "Forces (Hartrees/Bohr)" not in line:
                continue
            block = self._read_force_block(i)
            if block and last_input in self._input_frames:
                turned = self._torque_free(
                    self._input_frames[last_input], [g["vector"] for g in block]
                )
                if turned is None:
                    logging.warning(
                        "Gaussian: the forces at line %d have a net torque; "
                        "their orientation is approximate",
                        i + 1,
                    )
                else:
                    for g, vector in zip(block, turned):
                        g["vector"] = [float(v) for v in vector]
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
            input_coords = self._input_geometry(final["geometry_line"])[1]
        if input_coords and self.data.get("coords"):
            rotation = self._rotation_between(input_coords, self.data["coords"])
            if rotation is not None:
                for g in grads:
                    g["vector"] = [float(v) for v in rotation @ np.array(g["vector"])]
        self.data["gradients"] = grads

    def _input_geometry(self, header_idx):
        """(atoms, coords) of a table, in the input orientation where known."""
        atoms, coords = self._read_geometry(header_idx)
        frames = getattr(self, "_input_frames", {})
        return atoms, frames.get(header_idx, coords)

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
        starts = [0] + [i + 1 for i, l in enumerate(self.lines) if _ends_job(l)]
        starts = sorted(set(s for s in starts if s < len(self.lines)))
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
        ranges = self._job_ranges()
        irc = [r for r in ranges if any(_IRC_BANNER in l for l in self.lines[r[0] : r[1]])]
        if irc:
            # The reaction path is the result; the TS search before it is not.
            for start, end in irc:
                self._read_irc(start, end)
            return
        for start, end in ranges:
            job = self.lines[start:end]
            if any(_STEP.match(line) for line in job) and self._optimizes(job):
                self._read_optimization(start, end)
            elif any(line.strip() == _RIGID_SCAN for line in job):
                self._read_rigid_scan(start, end)

    @staticmethod
    def _optimizes(job):
        """False for a job step whose route asks for no optimization.

        The Freq step of Opt Freq checks the stationary point with one more
        "Step number" that is not a cycle. A job without a route is kept.
        """
        route = _route_text(job)
        return route is None or bool(_OPT_KEYWORD.search(route))

    def _read_rigid_scan(self, start, end):
        """Append one point per single-point energy of a rigid (Scan) job."""
        self.data["is_scan"] = True
        grads_at = {g["line"]: g["grads"] for g in self.data["all_gradients"]}
        geometry, points = None, []
        for i in range(start, end):
            line = self.lines[i]
            header = _GEOMETRY_HEADER.match(line)
            if header:
                if header.group(1) != "Standard" or geometry is None:
                    geometry = self._input_geometry(i)
                continue
            if i in grads_at:
                if points:
                    points[-1]["gradients"] = grads_at[i]
                continue
            scf = _SCF_DONE.match(line)
            if scf and geometry is not None:
                atoms, coords = geometry
                points.append(
                    {
                        "type": "scan_step",
                        "scan_step_id": len(points) + 1,
                        "step": 1,
                        "energy": float(scf.group(1)),
                        "atoms": atoms,
                        "coords": coords,
                        "convergence": {},
                        "gradients": [],
                    }
                )
                geometry = None
        values = self._rigid_scan_values(start, end)
        if len(values) == len(points):
            for point, value in zip(points, values):
                point["scan_coord"] = value
        self.data["scan_steps"].extend(points)

    def _rigid_scan_values(self, start, end):
        """First scanned variable of each point in the end-of-scan summary."""
        for i in range(start, end):
            if "Summary of the potential surface scan:" not in self.lines[i]:
                continue
            values = []
            for line in self.lines[i + 3 : end]:
                parts = line.split()
                if len(parts) < 3 or not parts[0].isdigit():
                    break
                try:
                    values.append(float(parts[1]))
                except ValueError:
                    break
            return values
        return []

    def _irc_structure(self, title_idx):
        """(atoms, coords) of the CURRENT STRUCTURE table below title_idx."""
        atoms, coords = [], []
        dashes = 0
        for line in self.lines[title_idx + 1 : title_idx + 400]:
            if line.strip().startswith("---"):
                dashes += 1
                if dashes == 3:
                    break
                continue
            parts = line.split()
            if dashes != 2 or len(parts) != 5:
                continue
            try:
                atoms.append(symbol_for(int(parts[1])))
                coords.append([float(v) for v in parts[2:5]])
            except ValueError:
                return [], []
        return atoms, coords

    def _checkpoint_structure(self, start, end):
        """(atoms, coords) of the "Sym,0,x,y,z" rows read from the checkpoint."""
        for i in range(start, end):
            line = self.lines[i]
            if "found in file" not in line and "Structure from the checkpoint" not in line:
                continue
            atoms, coords = [], []
            for row in self.lines[i + 1 : end]:
                parts = row.strip().split(",")
                try:
                    if len(parts) not in (4, 5):
                        raise ValueError(row)
                    xyz = [float(v) for v in parts[-3:]]
                except ValueError:
                    if atoms:
                        break
                    continue
                atoms.append(re.sub(r"[^A-Za-z].*$", "", parts[0]))
                coords.append(xyz)
            if atoms:
                return atoms, coords
        return None

    def _read_irc(self, start, end):
        """Append the reaction path, ordered from the reverse end to the forward end."""
        grads_at = {g["line"]: g["grads"] for g in self.data["all_gradients"]}
        energy, grads, geometry, direction = None, [], None, 1.0
        ts, path = None, []
        for i in range(start, end):
            line = self.lines[i]
            if _GEOMETRY_HEADER.match(line):
                geometry = self._input_geometry(i)
                continue
            scf = _SCF_DONE.match(line)
            if scf:
                energy = float(scf.group(1))
                continue
            if i in grads_at:
                grads = grads_at[i]
                continue
            match = _IRC_DIRECTION.match(line)
            if match:
                direction = -1.0 if match.group(1) == "REVERSE" else 1.0
                continue
            match = _IRC_POINT.match(line)
            if not match:
                continue
            if int(match.group(1)) == 0:
                if ts is None:
                    ts = self._irc_ts(start, i, geometry)
                continue
            if i + 1 >= end or "CURRENT STRUCTURE" not in self.lines[i + 1]:
                continue
            atoms, coords = self._irc_structure(i + 1)
            net = None
            for row in self.lines[i + 1 : min(i + 400, end)]:
                found = _IRC_NET.search(row)
                if found:
                    net = float(found.group(1))
                    break
            if atoms and energy is not None and net is not None:
                path.append((direction * net, energy, atoms, coords, grads))
        if ts is not None:
            if path and np is not None and len(ts[3]) == len(path[0][3]):
                # The checkpoint frame may be shifted from the path's.
                nearest = min(path, key=lambda p: abs(p[0]))
                ts = (0.0, ts[1], ts[2], self._best_fit(ts[3], nearest[3]), [])
            path.append(ts)
        path.sort(key=lambda p: p[0])
        for n, (coord, e, atoms, coords, g) in enumerate(path):
            self.data["scan_steps"].append(
                {
                    "type": "scan_step",
                    "scan_step_id": n + 1,
                    "step": 1,
                    "energy": e,
                    "atoms": atoms,
                    "coords": coords,
                    "convergence": {},
                    "gradients": g,
                    "scan_coord": coord,
                }
            )

    def _irc_ts(self, start, point_line, geometry):
        """(0.0, energy, atoms, coords, []) of the transition state, or None."""
        chk_energy = None
        for line in self.lines[start:point_line]:
            found = _IRC_CHK_ENERGY.match(line)
            if found:
                chk_energy = float(found.group(1))
        energy = None
        for line in reversed(self.lines[:point_line]):
            scf = _SCF_DONE.match(line)
            if scf:
                energy = float(scf.group(1))
                break
        if chk_energy is not None and (energy is None or abs(energy - chk_energy) > 1e-6):
            energy = chk_energy
        if geometry is None:
            geometry = self._checkpoint_structure(start, point_line)
        if geometry is None or energy is None:
            return None
        return (0.0, energy, geometry[0], geometry[1], [])

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
                    current["geometry"] = self._input_geometry(i)
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
