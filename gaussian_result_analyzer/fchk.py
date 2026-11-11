"""Gaussian formatted checkpoint files.

A .fchk holds what a log only prints on request: the basis set and every orbital,
at full precision. It holds no SCF trace, thermochemistry, NMR, NBO or Hirshfeld
results, so a log beside it still adds those. Lengths in the file are in bohr.
"""

import re

BOHR_TO_ANGSTROM = 0.52917721067
AU_TO_DEBYE = 2.541746473
HARTREE_TO_EV = 27.211386245988
HARTREE_TO_CM = 219474.6313632
EV_NM = 1239.8419843

#: The extensions formchk writes: .fchk elsewhere, .fch by default on Windows.
FCHK_EXTENSIONS = (".fchk", ".fch")

#: A log and a checkpoint belong together only within these limits.
MAX_COORD_DIFFERENCE_ANGSTROM = 2e-3
MAX_ENERGY_DIFFERENCE_HARTREE = 1e-5

_HEADER = re.compile(r"^(\S.{0,39}?)\s+([IRCL])(\s+N=)?\s+(\S+)\s*$")

#: Sections read; the rest of a large file is skipped without converting it.
_WANTED = {
    "Number of atoms",
    "Route",
    "Charge",
    "Multiplicity",
    "Number of alpha electrons",
    "Number of beta electrons",
    "Atomic numbers",
    "Current cartesian coordinates",
    "Shell types",
    "Number of primitives per shell",
    "Shell to atom map",
    "Primitive exponents",
    "Contraction coefficients",
    "P(S=P) Contraction coefficients",
    "SCF Energy",
    "Total Energy",
    "MP2 Energy",
    "Alpha Orbital Energies",
    "Beta Orbital Energies",
    "Alpha MO coefficients",
    "Beta MO coefficients",
    "Mulliken Charges",
    "NPA Charges",
    "ESP Charges",
    "Nuclear derivative order",
    "Cartesian Gradient",
    "Dipole Moment",
    "S**2",
    "Number of Normal Modes",
    "Vib-E2",
    "Vib-Modes",
    "Optimization Num results per geometry",
    "Optimization Num geometry variables",
    "Optimization Number of geometries",
    "ETran scalars",
    "ETran state values",
    "Gaussian Version",
}
_OPT_POINT = re.compile(r"^Opt point\s+(\d+) (Results|Geometries|Gradient)")

_SYMBOLS = (
    "X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu "
    "Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba "
    "La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi "
    "Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr"
).split()

_CARTESIAN_LABELS = {
    2: ("XX", "YY", "ZZ", "XY", "XZ", "YZ"),
    3: ("XXX", "YYY", "ZZZ", "XYY", "XXY", "XXZ", "XZZ", "YZZ", "YYZ", "XYZ"),
}
_SHELL_LETTER = "SPDFGHI"


def looks_like_fchk(text):
    """Whether the start of *text* has the layout of a formatted checkpoint."""
    lines = text.splitlines()[:6]
    return len(lines) >= 3 and any(
        line.startswith("Number of atoms") and _HEADER.match(line) for line in lines[2:]
    )


def read_sections(text):
    """The wanted sections as {name: scalar or list}; numbers are converted."""
    if not looks_like_fchk(text):
        raise ValueError("Not a Gaussian formatted checkpoint file.")
    lines = text.splitlines()
    sections = {}
    i = 2
    while i < len(lines):
        match = _HEADER.match(lines[i])
        i += 1
        if not match:
            continue
        name, kind, is_array, value = match.groups()
        start = i
        if is_array and kind == "C":
            # Text rows hold five 12-character words and may start with anything.
            i += (int(value) + 4) // 5
        elif is_array:
            while i < len(lines) and lines[i][:1] == " ":
                i += 1
        if name not in _WANTED and not _OPT_POINT.match(name):
            continue
        convert = {"I": int, "R": float}.get(kind, str)
        try:
            if is_array:
                if kind == "C":
                    sections[name] = "".join(lines[start:i]).strip()
                    continue
                tokens = " ".join(lines[start:i]).split()
                if len(tokens) != int(value):
                    raise ValueError(f"expected {value} values, found {len(tokens)}")
                sections[name] = [convert(token) for token in tokens]
            else:
                sections[name] = convert(value)
        except ValueError as exc:
            raise ValueError(f"Damaged section '{name}': {exc}") from exc
    if "Atomic numbers" not in sections or (
        "Current cartesian coordinates" not in sections
    ):
        raise ValueError("The checkpoint has no molecular structure.")
    return sections


def _symbol(number):
    return _SYMBOLS[number] if 0 <= number < len(_SYMBOLS) else "X"


def _triples(values, factor=1.0):
    return [
        [values[k] * factor, values[k + 1] * factor, values[k + 2] * factor]
        for k in range(0, len(values) - 2, 3)
    ]


def structure(sections):
    """(element symbols, coordinates in Angstrom)."""
    atoms = [_symbol(number) for number in sections["Atomic numbers"]]
    if len(atoms) != sections.get("Number of atoms", len(atoms)) or len(
        sections["Current cartesian coordinates"]
    ) != 3 * len(atoms):
        raise ValueError("The checkpoint coordinates do not match the number of atoms.")
    return atoms, _triples(
        sections["Current cartesian coordinates"], BOHR_TO_ANGSTROM
    )


def basis_shells(sections, coords):
    """Contracted shells, with an SP shell split into an S and a P."""
    types = sections.get("Shell types")
    counts = sections.get("Number of primitives per shell")
    atoms = sections.get("Shell to atom map")
    exponents = sections.get("Primitive exponents")
    coefficients = sections.get("Contraction coefficients")
    if not (types and counts and atoms and exponents and coefficients):
        return []
    if not (len(types) == len(counts) == len(atoms)) or any(
        count < 1 for count in counts
    ):
        raise ValueError("The checkpoint shell arrays have inconsistent lengths.")
    if any(atom < 1 or atom > len(coords) for atom in atoms):
        raise ValueError("The checkpoint shell map refers to a missing atom.")
    if sum(counts) != len(exponents) or len(exponents) != len(coefficients):
        raise ValueError("The checkpoint primitive arrays have inconsistent lengths.")
    sp_coefficients = sections.get("P(S=P) Contraction coefficients")
    if -1 in types and (not sp_coefficients or len(sp_coefficients) != len(exponents)):
        raise ValueError("An SP shell has incomplete P coefficients.")
    if any(abs(shell_type) >= len(_SHELL_LETTER) for shell_type in types):
        raise ValueError("The checkpoint contains an unsupported shell type.")
    shells = []
    start = 0
    for shell_type, count, atom in zip(types, counts, atoms):
        end = start + count
        atom_idx = atom - 1
        common = {
            "atom_idx": atom_idx,
            "origin": coords[atom_idx],
            "exps": exponents[start:end],
        }
        if shell_type == -1:
            if not sp_coefficients:
                raise ValueError("An SP shell has no P coefficients.")
            shells.append(
                dict(common, l=0, coeffs=coefficients[start:end], cartesian=False)
            )
            shells.append(
                dict(common, l=1, coeffs=sp_coefficients[start:end], cartesian=False)
            )
        else:
            shells.append(
                dict(
                    common,
                    l=abs(shell_type),
                    coeffs=coefficients[start:end],
                    cartesian=shell_type >= 2,
                )
            )
        start = end
    return shells


def function_labels(sections, atoms):
    """(atom index, element, function label) for each AO, in Gaussian's order."""
    labels = []
    numbers = {}
    for shell_type, atom in zip(
        sections.get("Shell types", []), sections.get("Shell to atom map", [])
    ):
        atom_idx = atom - 1
        n = numbers[atom_idx] = numbers.get(atom_idx, 0) + 1
        if shell_type == 0:
            names = ("S",)
        elif shell_type == 1:
            names = ("PX", "PY", "PZ")
        elif shell_type == -1:
            names = ("S", "PX", "PY", "PZ")
        elif shell_type in _CARTESIAN_LABELS:
            names = _CARTESIAN_LABELS[shell_type]
        elif shell_type < 0:
            letter = _SHELL_LETTER[-shell_type]
            names = [f"{letter} 0"] + [
                f"{letter}{sign}{m}"
                for m in range(1, -shell_type + 1)
                for sign in "+-"
            ]
        else:
            count = (shell_type + 1) * (shell_type + 2) // 2
            letter = _SHELL_LETTER[shell_type]
            names = [f"{letter}{k + 1}" for k in range(count)]
        labels.extend((atom_idx, atoms[atom_idx], f"{n}{name}") for name in names)
    return labels


def molecular_orbitals(sections, atoms):
    """data["mo_coeffs"] from the alpha and, when present, beta tables."""
    n_alpha = sections.get("Number of alpha electrons", 0)
    n_beta = sections.get("Number of beta electrons", 0)
    labels = function_labels(sections, atoms)
    n_basis = len(labels)
    unrestricted = "Beta MO coefficients" in sections
    tables = (
        (
            ("alpha", "Alpha", lambda i: float(i < n_alpha)),
            ("beta", "Beta", lambda i: float(i < n_beta)),
        )
        if unrestricted
        else (
            (
                "restricted",
                "Alpha",
                lambda i: float(i < n_alpha) + float(i < n_beta),
            ),
        )
    )
    orbitals = {}
    for spin, name, occupation in tables:
        energies = sections.get(f"{name} Orbital Energies")
        values = sections.get(f"{name} MO coefficients")
        if not energies or not values or not n_basis:
            return {}
        if len(values) != len(energies) * n_basis:
            raise ValueError(
                f"{len(values)} {name} MO coefficients do not fill "
                f"{len(energies)} orbitals of {n_basis} functions."
            )
        for index, energy in enumerate(energies):
            row = values[index * n_basis : (index + 1) * n_basis]
            orbitals[f"{index}_{spin}"] = {
                "coeffs": [
                    {"atom_idx": atom, "sym": symbol, "orb": orb, "coeff": coeff}
                    for (atom, symbol, orb), coeff in zip(labels, row)
                ],
                "spin": spin,
                "energy": energy,
                "occ": occupation(index),
                "id": index,
            }
    return orbitals


def normal_modes(sections, n_atoms):
    """data["frequencies"]: Vib-E2 holds 14 quantities for each mode, in rows."""
    n_modes = sections.get("Number of Normal Modes")
    e2 = sections.get("Vib-E2")
    modes = sections.get("Vib-Modes")
    if not n_modes or not e2 or not modes:
        return []
    if len(modes) != n_modes * 3 * n_atoms or len(e2) < 5 * n_modes:
        raise ValueError("The normal modes do not match the number of atoms.")
    return [
        {
            "freq": e2[k],
            "reduced_mass": e2[n_modes + k],
            "force_constant": e2[2 * n_modes + k],
            "ir": e2[3 * n_modes + k],
            "raman": e2[4 * n_modes + k],
            "vector": _triples(modes[3 * n_atoms * k : 3 * n_atoms * (k + 1)]),
        }
        for k in range(n_modes)
    ]


def gradients(sections, atoms):
    """data["gradients"], only when the job computed them: a single point has zeros."""
    values = sections.get("Cartesian Gradient")
    if not values or sections.get("Nuclear derivative order", 0) < 1:
        return []
    return [
        {"atom_idx": k, "atom_sym": atoms[k], "vector": vector}
        for k, vector in enumerate(_triples(values))
    ]


def dipole(sections):
    """data["dipole"] from the vector in atomic units."""
    vector = sections.get("Dipole Moment")
    if not vector or len(vector) != 3:
        return None
    magnitude = sum(v * v for v in vector) ** 0.5
    debye = [v * AU_TO_DEBYE for v in vector]
    return {
        "vector_au": vector,
        "vector_debye": debye,
        "magnitude_au": magnitude,
        "magnitude_debye": magnitude * AU_TO_DEBYE,
        "vector": debye,
        "magnitude": magnitude * AU_TO_DEBYE,
    }


def _esp_label(route):
    """The name the log gives these fitted charges, read from the route."""
    upper = route.upper()
    if "CHELPG" in upper:
        return "CHELPG"
    if re.search(r"\bMK\b", upper):
        return "Merz-Kollman"
    return "ESP"


def atomic_charges(sections, atoms):
    """data["charges"]; the checkpoint has no APT, Hirshfeld or CM5 charges."""
    charges = {}
    for label, name in (
        ("Mulliken", "Mulliken Charges"),
        ("NBO", "NPA Charges"),
        (_esp_label(sections.get("Route", "")), "ESP Charges"),
    ):
        values = sections.get(name)
        if values and len(values) == len(atoms):
            charges[label] = [
                {"atom_idx": k, "atom_sym": atoms[k], "charge": value}
                for k, value in enumerate(values)
            ]
    return charges


def optimization_steps(sections, atoms):
    """data["scan_steps"]: more than one point means a scan."""
    per_point = sections.get("Optimization Number of geometries")
    n_results = sections.get("Optimization Num results per geometry")
    n_variables = sections.get("Optimization Num geometry variables")
    if not per_point or not n_results or not n_variables or sum(per_point) < 2:
        return [], False
    scan = len(per_point) > 1
    steps = []
    for point, count in enumerate(per_point, 1):
        prefix = f"Opt point {point:7d}"
        results = sections.get(f"{prefix} Results for each geome", [])
        geometries = sections.get(f"{prefix} Geometries", [])
        grads = sections.get(f"{prefix} Gradient at each geome", [])
        for k in range(count):
            geometry = geometries[k * n_variables : (k + 1) * n_variables]
            grad = grads[k * n_variables : (k + 1) * n_variables]
            if len(geometry) != n_variables or 2 * k >= len(results):
                raise ValueError(f"Optimization point {point} is incomplete.")
            steps.append(
                {
                    "type": "scan_step" if scan else "opt_cycle",
                    "scan_step_id": point,
                    "step": len(steps) + 1,
                    "energy": results[n_results * k],
                    "atoms": list(atoms),
                    "coords": _triples(geometry, BOHR_TO_ANGSTROM),
                    "convergence": {},
                    "gradients": [
                        {"atom_idx": a, "atom_sym": atoms[a], "vector": v}
                        for a, v in enumerate(_triples(grad))
                    ],
                }
            )
    return steps, scan


def excited_states(sections):
    """data["tddft"]: vertical energies and length-gauge oscillator strengths."""
    scalars = sections.get("ETran scalars")
    values = sections.get("ETran state values")
    reference = sections.get("SCF Energy")
    if not scalars or not values or reference is None or len(scalars) < 2:
        return []
    n_states, stride = scalars[0], scalars[1]
    if len(values) < n_states * stride or stride < 4:
        return []
    states = []
    for k in range(n_states):
        row = values[k * stride : (k + 1) * stride]
        excitation = row[0] - reference
        ev = excitation * HARTREE_TO_EV
        strength = 2.0 / 3.0 * excitation * sum(v * v for v in row[1:4])
        states.append(
            {
                "state": k + 1,
                "label": "",
                "energy_ev": ev,
                "energy_nm": EV_NM / ev if ev > 0 else 0.0,
                "energy_cm": excitation * HARTREE_TO_CM,
                "osc": strength,
                "osc_len": strength,
                "osc_vel": 0.0,
                "rotatory_strength": 0.0,
                "rot_len": 0.0,
                "rot_vel": 0.0,
                "transitions": [],
            }
        )
    return states


def why_it_differs(log_data, sections):
    """Why a checkpoint does not belong to a parsed log, or "" when it does."""
    atoms, coords = structure(sections)
    if atoms != log_data.get("atoms"):
        return "its atoms differ from the output file's"
    for key, name in (("charge", "Charge"), ("mult", "Multiplicity")):
        if name in sections and sections[name] != log_data.get(key):
            return f"its {name.lower()} differs from the output file's"
    log_coords = log_data.get("coords") or []
    if len(coords) != len(log_coords) or any(len(row) != 3 for row in log_coords):
        return "its geometry or orientation differs from the output file's"
    largest = max(
        (
            abs(a - b)
            for pair in zip(coords, log_coords)
            for a, b in zip(*pair)
        ),
        default=float("inf"),
    )
    if largest > MAX_COORD_DIFFERENCE_ANGSTROM:
        return "its geometry or orientation differs from the output file's"
    log_energy = log_data.get("scf_energy")
    fchk_energy = sections.get("SCF Energy")
    if (
        log_energy is not None
        and fchk_energy is not None
        and abs(log_energy - fchk_energy) > MAX_ENERGY_DIFFERENCE_HARTREE
    ):
        return "its SCF energy differs from the output file's"
    return ""
