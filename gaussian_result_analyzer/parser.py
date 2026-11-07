"""GaussianParser: assembles the parsing mixins and drives parse_all's step order."""

from .parser_structure import _StructureParsingMixin
from . import fchk
from .parser_electronic import _ElectronicParsingMixin
from .parser_properties import _PropertyParsingMixin, AU_TO_DEBYE
from .parser_spectra import _SpectraParsingMixin, IMAGINARY_FREQ_THRESHOLD

__all__ = [
    "ParseCancelled",
    "GaussianParser",
    "AU_TO_DEBYE",
    "IMAGINARY_FREQ_THRESHOLD",
]


class ParseCancelled(Exception):
    """Raised by a progress callback to abort ``parse_all`` between steps."""


class GaussianParser(
    _StructureParsingMixin,
    _ElectronicParsingMixin,
    _PropertyParsingMixin,
    _SpectraParsingMixin,
):
    """Parser for Gaussian output files (.log, .out)"""

    def __init__(self):
        self.filename = ""
        self.raw_content = ""
        self.lines = []
        self._reset_data()

    def _reset_data(self):
        """Discard every result derived from the previous input."""
        self.data = {
            "scf_traces": [],  # List of SCF energy values per iteration
            "converged": False,
            "scf_energy": None,
            "atoms": [],
            "coords": [],
            "charge": 0,
            "mult": 1,
            "frequencies": [],  # List of dicts: {freq, ir, raman, vector}
            "excitation_energies": [],  # TDDFT
            "dipole": None,
            "dipoles": None,
            "nmr_shielding": [],
            "nmr_couplings": [],
            "charges": {},  # Type -> List
            "version": None,
            "scan_steps": [],
            "all_steps": [],
            "termination_status": "Running",
            "fchk_path": "",
            "fchk_note": "",
        }

    #: (method name, label) in execution order. Gradients must be parsed
    #: before the trajectory so they can be linked to it.
    PARSE_STEPS = (
        ("parse_basic", "Reading job information"),
        ("parse_gradients", "Reading forces"),
        ("parse_trajectory", "Reading optimization steps"),
        ("parse_frequencies", "Reading vibrational frequencies"),
        ("parse_thermal", "Reading thermochemistry"),
        ("parse_mo_coeffs", "Reading MO coefficients"),
        ("parse_orbital_energies", "Reading orbital energies"),
        ("parse_charges", "Reading atomic charges"),
        ("parse_nbo_orbitals", "Reading NBO orbitals"),
        ("parse_nbo_perturbation", "Reading NBO perturbation analysis"),
        ("parse_dipole", "Reading dipole moment"),
        ("parse_spin_contamination", "Reading spin contamination"),
        ("parse_dispersion", "Reading dispersion correction"),
        ("parse_energy_components", "Reading energy components"),
        ("parse_tddft", "Reading excited states"),
        ("parse_nmr", "Reading NMR data"),
        ("parse_basis_set", "Reading basis set"),
        ("parse_scf_trace", "Reading SCF convergence"),
        ("parse_scan_results_table", "Reading scan results"),
    )

    def load_from_memory(self, content, filename="", progress=None):
        """Store the Gaussian output text and run parse_all on it."""
        self._reset_data()
        self.filename = filename
        self.raw_content = content
        self.lines = content.splitlines()
        self.parse_all(progress=progress)

    def load_fchk(self, content, filename="", progress=None):
        """Fill the data from a formatted checkpoint, with no output file."""
        self._reset_data()
        self.filename = filename
        self.raw_content = ""
        self.lines = []
        if progress is not None:
            progress(0, 1, "Reading formatted checkpoint")
        sections = fchk.read_sections(content)
        atoms, coords = fchk.structure(sections)
        data = self.data
        data["atoms"] = atoms
        data["coords"] = coords
        data["charge"] = sections.get("Charge", 0)
        data["mult"] = sections.get("Multiplicity", 1)
        data["scf_energy"] = sections.get("SCF Energy", sections.get("Total Energy"))
        data["version"] = sections.get("Gaussian Version")
        data["termination_status"] = "Formatted checkpoint"
        data["charges"] = fchk.atomic_charges(sections, atoms)
        data["dipole"] = data["dipoles"] = fchk.dipole(sections)
        data["gradients"] = fchk.gradients(sections, atoms)
        data["frequencies"] = fchk.normal_modes(sections, len(atoms))
        data["scan_steps"], data["is_scan"] = fchk.optimization_steps(sections, atoms)
        data["tddft"] = fchk.excited_states(sections)
        s2 = sections.get("S**2")
        if s2 is not None:
            spin = (data["mult"] - 1) / 2
            ideal = spin * (spin + 1)
            data["spin_s2"] = {
                "actual": s2,
                "ideal": ideal,
                "contamination": s2 - ideal,
            }
        if "MP2 Energy" in sections:
            data["energy_components"] = [
                {
                    "label": "MP2 Total",
                    "value": sections["MP2 Energy"],
                    "dimensionless": False,
                }
            ]
        self._take_orbitals(sections)
        if progress is not None:
            progress(1, 1, "Finished")

    def merge_fchk(self, content):
        """Take the orbitals from a checkpoint of this job; say why not otherwise."""
        sections = fchk.read_sections(content)
        reason = fchk.why_it_differs(self.data, sections)
        if not reason:
            self._take_orbitals(sections)
        return reason

    def _take_orbitals(self, sections):
        """Replace orbitals and basis with the checkpoint's, when it has them."""
        atoms, coords = fchk.structure(sections)
        orbitals = fchk.molecular_orbitals(sections, atoms)
        if not orbitals:
            return
        types = sections["Shell types"]
        self.data["mo_coeffs"] = orbitals
        self.data["cartesian_d"] = 2 in types
        self.data["cartesian_f"] = 3 in types
        self.data["basis_set_shells"] = fchk.basis_shells(sections, coords)
        self.parse_orbital_energies()

    def parse_all(self, progress=None):
        """Run every parse step.

        *progress* is called as ``progress(done, total, label)`` before each
        step and once when finished; it may raise ``ParseCancelled`` to abort.
        """
        total = len(self.PARSE_STEPS)
        for done, (method, label) in enumerate(self.PARSE_STEPS):
            if progress is not None:
                progress(done, total, label)
            getattr(self, method)()
        if progress is not None:
            progress(total, total, "Finished")
