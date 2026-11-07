"""Formatted checkpoints against the log of the same Gaussian 16 job.

Every pair in tests/sample_outputs comes from one job, so the two files must agree
to the precision the log prints. The orbital test is the independent one: orbitals
read from a checkpoint alone, evaluated on a real-space grid, must be orthonormal.
"""

import importlib
import os
import sys
import types
import unittest
from unittest.mock import patch

try:
    import numpy as np
except ImportError:
    np = None

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gui_harness  # noqa: E402

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_SAMPLES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sample_outputs")
_PKG = "_gaussian_fchk_only"


def _module(name):
    if _PKG not in sys.modules:
        package = types.ModuleType(_PKG)
        package.__path__ = [os.path.join(_ROOT, "gaussian_result_analyzer")]
        package.PLUGIN_VERSION = "test"
        sys.modules[_PKG] = package
    return importlib.import_module(f"{_PKG}.{name}")


def _text(name):
    with open(os.path.join(_SAMPLES, name), encoding="utf-8") as handle:
        return handle.read()


def _log(name):
    parser = _module("parser").GaussianParser()
    parser.load_from_memory(_text(name + ".out"), name + ".out")
    return parser


def _fchk(name):
    parser = _module("parser").GaussianParser()
    parser.load_fchk(_text(name + ".fchk"), name + ".fchk")
    return parser


class TestAgainstTheLog(unittest.TestCase):
    def test_structure_charge_and_energy(self):
        for name in ("water_fk_mo", "o2_fk_uhf", "benzene_fk_opt_freq"):
            with self.subTest(job=name):
                log, fchk = _log(name).data, _fchk(name).data
                self.assertEqual(fchk["atoms"], log["atoms"])
                self.assertEqual(
                    (fchk["charge"], fchk["mult"]), (log["charge"], log["mult"])
                )
                self.assertAlmostEqual(fchk["scf_energy"], log["scf_energy"], places=8)
                for row, other in zip(fchk["coords"], log["coords"]):
                    for a, b in zip(row, other):
                        self.assertAlmostEqual(a, b, places=5)

    def test_orbitals_are_the_logs_at_higher_precision(self):
        for name in (
            "water_fk_mo",
            "h2co_fk_5d7f",
            "h2co_fk_cart",
            "o2_fk_uhf",
            "o2_fk_rohf",
        ):
            with self.subTest(job=name):
                log, fchk = _log(name).data, _fchk(name).data
                self.assertEqual(set(fchk["mo_coeffs"]), set(log["mo_coeffs"]))
                for key, orbital in log["mo_coeffs"].items():
                    other = fchk["mo_coeffs"][key]
                    self.assertAlmostEqual(other["energy"], orbital["energy"], delta=6e-6)
                    self.assertEqual(other["occ"], orbital["occ"])
                    self.assertEqual(len(other["coeffs"]), len(orbital["coeffs"]))
                    for a, b in zip(orbital["coeffs"], other["coeffs"]):
                        self.assertAlmostEqual(a["coeff"], b["coeff"], delta=6e-6)
                        self.assertEqual(
                            (a["atom_idx"], a["sym"], a["orb"].strip()),
                            (b["atom_idx"], b["sym"], b["orb"].strip()),
                        )

    def test_shells_are_the_logs(self):
        for name in ("water_fk_mo", "h2co_fk_5d7f", "h2co_fk_cart", "o2_fk_uhf"):
            with self.subTest(job=name):
                log, fchk = _log(name).data, _fchk(name).data
                self.assertEqual(
                    (fchk["cartesian_d"], fchk["cartesian_f"]),
                    (log["cartesian_d"], log["cartesian_f"]),
                )
                shells = log["basis_set_shells"]
                self.assertEqual(len(fchk["basis_set_shells"]), len(shells))
                for a, b in zip(shells, fchk["basis_set_shells"]):
                    self.assertEqual(
                        (a["l"], a["atom_idx"], a["cartesian"]),
                        (b["l"], b["atom_idx"], b["cartesian"]),
                    )
                    for x, y in zip(a["exps"], b["exps"]):
                        self.assertAlmostEqual(x / y, 1.0, places=6)
                    for x, y in zip(a["coeffs"], b["coeffs"]):
                        self.assertAlmostEqual(x, y, places=6)

    def test_normal_modes(self):
        for name in ("water_fk_raman", "benzene_fk_opt_freq"):
            with self.subTest(job=name):
                log, fchk = _log(name).data, _fchk(name).data
                self.assertEqual(len(fchk["frequencies"]), len(log["frequencies"]))
                for a, b in zip(log["frequencies"], fchk["frequencies"]):
                    for key in ("freq", "ir", "raman", "reduced_mass", "force_constant"):
                        self.assertAlmostEqual(a[key], b[key], delta=1.1e-3)
                    for p, q in zip(a["vector"], b["vector"]):
                        for x, y in zip(p, q):
                            self.assertAlmostEqual(x, y, delta=5.1e-3)

    def test_gradients(self):
        for name in ("water_fk_force", "water_fk_opt", "benzene_fk_opt_freq"):
            with self.subTest(job=name):
                log, fchk = _log(name).data, _fchk(name).data
                self.assertEqual(len(fchk["gradients"]), len(log["gradients"]))
                for a, b in zip(log["gradients"], fchk["gradients"]):
                    for x, y in zip(a["vector"], b["vector"]):
                        self.assertAlmostEqual(x, y, places=7)

    def test_a_single_point_has_no_gradients_though_the_section_is_there(self):
        self.assertIn("Cartesian Gradient", _text("water_fk_mo.fchk"))
        self.assertEqual(_fchk("water_fk_mo").data["gradients"], [])

    def test_dipole(self):
        for name in ("water_fk_mo", "h2co_fk_5d7f"):
            with self.subTest(job=name):
                log, fchk = _log(name).data["dipole"], _fchk(name).data["dipole"]
                self.assertAlmostEqual(
                    log["magnitude_debye"], fchk["magnitude_debye"], places=3
                )
                for x, y in zip(log["vector_debye"], fchk["vector_debye"]):
                    self.assertAlmostEqual(x, y, places=3)

    def test_charges(self):
        for name, kinds in (
            ("water_fk_mo", ("Mulliken",)),
            ("water_fk_raman", ("Mulliken", "NBO")),
            ("water_fk_chelpg", ("Mulliken", "CHELPG")),
        ):
            log, fchk = _log(name).data["charges"], _fchk(name).data["charges"]
            self.assertEqual(set(fchk), set(kinds))
            for kind in kinds:
                with self.subTest(job=name, kind=kind):
                    for a, b in zip(log[kind], fchk[kind]):
                        self.assertEqual(a["atom_sym"], b["atom_sym"])
                        self.assertAlmostEqual(a["charge"], b["charge"], places=3)

    def test_the_log_keeps_the_charges_of_the_last_geometry(self):
        """Pop=NBO in an Opt prints the table at each end; the file was read first."""
        npa = _log("water_fk_raman").data["charges"]["NBO"]
        self.assertAlmostEqual(npa[0]["charge"], -0.95191)

    def test_excited_states(self):
        log, fchk = _log("h2co_fk_tddft").data, _fchk("h2co_fk_tddft").data
        self.assertEqual(len(fchk["tddft"]), len(log["tddft"]))
        for a, b in zip(log["tddft"], fchk["tddft"]):
            self.assertAlmostEqual(a["energy_ev"], b["energy_ev"], places=3)
            self.assertAlmostEqual(a["energy_nm"], b["energy_nm"], delta=0.01)
            self.assertAlmostEqual(a["osc"], b["osc"], delta=6e-5)

    def test_optimization_and_scan(self):
        for name, scan in (("water_fk_opt", False), ("ethane_fk_scan", True)):
            with self.subTest(job=name):
                log, fchk = _log(name).data, _fchk(name).data
                self.assertEqual(bool(log.get("is_scan")), scan)
                self.assertEqual(fchk["is_scan"], scan)
                self.assertEqual(len(fchk["scan_steps"]), len(log["scan_steps"]))
                for a, b in zip(log["scan_steps"], fchk["scan_steps"]):
                    self.assertEqual(a["type"], b["type"])
                    self.assertEqual(a["scan_step_id"], b["scan_step_id"])
                    self.assertAlmostEqual(a["energy"], b["energy"], places=6)

    def test_spin_and_post_hf_energy(self):
        s2 = _fchk("o2_fk_uhf").data["spin_s2"]
        self.assertAlmostEqual(s2["actual"], _log("o2_fk_uhf").data["spin_s2"]["actual"], 3)
        total = {
            row["label"]: row["value"] for row in _fchk("water_fk_mp2").data["energy_components"]
        }
        logged = {
            row["label"]: row["value"] for row in _log("water_fk_mp2").data["energy_components"]
        }
        self.assertAlmostEqual(total["MP2 Total"], logged["MP2 Total"], places=8)


class TestWhatOnlyTheCheckpointHolds(unittest.TestCase):
    def test_a_plain_job_prints_no_orbitals_but_its_checkpoint_has_them(self):
        self.assertEqual(_log("water_fk_plain").data.get("mo_coeffs"), {})
        self.assertEqual(len(_fchk("water_fk_plain").data["mo_coeffs"]), 19)

    def test_it_holds_no_scf_trace_or_thermochemistry(self):
        data = _fchk("benzene_fk_opt_freq").data
        self.assertEqual(data["scf_traces"], [])
        self.assertFalse(data.get("thermal"))
        self.assertEqual(data["nmr_shielding"], [])

    def test_a_single_geometry_makes_no_trajectory(self):
        self.assertEqual(_fchk("water_fk_mo").data["scan_steps"], [])


@unittest.skipIf(np is None, "numpy is required for the orbital integrals")
class TestOrbitalsFromACheckpointAreOrthonormal(unittest.TestCase):
    """Fails if an order, a normalization or a Cartesian/pure assignment is wrong."""

    def _assert_orthonormal(self, name, spin, count):
        from test_gaussian_parser import _grid

        data = _fchk(name).data
        engine_module = gui_harness.load_isolated("mo_engine")
        shells = [
            {
                "type": s["l"],
                "center": np.array(s["origin"]) / 0.529177249,
                "exps": np.array(s["exps"]),
                "coeffs": np.array(s["coeffs"]),
                "cartesian": bool(s["cartesian"]),
            }
            for s in data["basis_set_shells"]
        ]
        engine = engine_module.BasisSetEngine(shells)
        points, weights = _grid(np.array(data["coords"]) / 0.529177249)
        values = np.array(
            [
                engine.evaluate_mo_on_grid(
                    index,
                    points,
                    np.array(
                        [c["coeff"] for c in data["mo_coeffs"][f"{index}_{spin}"]["coeffs"]]
                    ),
                )
                for index in range(count)
            ]
        )
        overlap = (values * weights) @ values.T
        self.assertLess(np.abs(overlap - np.eye(count)).max(), 2e-2)

    def test_cartesian_d_with_sp_shells(self):
        self._assert_orthonormal("water_fk_plain", "restricted", 19)

    def test_pure_d_and_f(self):
        self._assert_orthonormal("h2co_fk_5d7f", "restricted", 30)

    def test_cartesian_d_and_f(self):
        self._assert_orthonormal("h2co_fk_cart", "restricted", 30)

    def test_pure_g(self):
        self._assert_orthonormal("water_fk_pvqz", "restricted", 30)

    def test_unrestricted(self):
        self._assert_orthonormal("o2_fk_uhf", "alpha", 30)
        self._assert_orthonormal("o2_fk_uhf", "beta", 30)

    def test_restricted_open_shell(self):
        self._assert_orthonormal("o2_fk_rohf", "restricted", 30)


class TestALogWithACheckpoint(unittest.TestCase):
    def test_charge_spin_and_incomplete_geometry_are_refused(self):
        for key, value, reason in (
            ("charge", 1, "charge"),
            ("mult", 3, "multiplicity"),
            ("coords", [[0.0, 0.0, 0.0]], "geometry"),
        ):
            with self.subTest(field=key):
                parser = _log("water_fk_plain")
                parser.data[key] = value
                self.assertIn(reason, parser.merge_fchk(_text("water_fk_plain.fchk")))
                self.assertEqual(parser.data["mo_coeffs"], {})

    def test_missing_orbitals_are_not_reported_as_a_success(self):
        parser = _log("water_fk_mo")
        original = parser.data["mo_coeffs"]
        sections = _module("fchk").read_sections(_text("water_fk_mo.fchk"))
        del sections["Alpha MO coefficients"]
        with patch.object(_module("fchk"), "read_sections", return_value=sections):
            self.assertIn("no complete", parser.merge_fchk("checkpoint"))
        self.assertIs(parser.data["mo_coeffs"], original)

    def test_damaged_basis_is_rejected_without_replacing_existing_orbitals(self):
        for name, value in (
            ("Shell to atom map", [99]),
            ("Shell to atom map", [0]),
            ("Number of primitives per shell", [0]),
            ("Primitive exponents", [1.0]),
            ("P(S=P) Contraction coefficients", [1.0]),
            ("Shell types", [99]),
        ):
            with self.subTest(section=name, value=value):
                parser = _log("water_fk_mo")
                original = parser.data["mo_coeffs"]
                sections = _module("fchk").read_sections(_text("water_fk_mo.fchk"))
                if name in ("Shell to atom map", "Number of primitives per shell", "Shell types"):
                    sections[name][0] = value[0]
                else:
                    sections[name] = value
                with patch.object(_module("fchk"), "read_sections", return_value=sections):
                    with self.assertRaises(ValueError):
                        parser.merge_fchk("checkpoint")
                self.assertIs(parser.data["mo_coeffs"], original)

    def test_the_orbitals_join_the_logs_results(self):
        parser = _log("water_fk_plain")
        self.assertEqual(parser.merge_fchk(_text("water_fk_plain.fchk")), "")
        self.assertEqual(len(parser.data["mo_coeffs"]), 19)
        self.assertEqual(len(parser.data["orbital_energies"]), 19)
        self.assertTrue(parser.data["basis_set_shells"])
        self.assertTrue(parser.data["scf_traces"])

    def test_another_molecule_is_refused(self):
        parser = _log("water_fk_plain")
        reason = parser.merge_fchk(_text("h2co_fk_5d7f.fchk"))
        self.assertIn("atoms", reason)
        self.assertEqual(parser.data["mo_coeffs"], {})

    def test_another_geometry_is_refused(self):
        parser = _log("water_fk_plain")
        reason = parser.merge_fchk(_text("water_fk_raman.fchk"))
        self.assertIn("geometry", reason)
        self.assertEqual(parser.data["mo_coeffs"], {})

    def test_another_method_at_the_same_geometry_is_refused(self):
        parser = _log("water_fk_plain")
        reason = parser.merge_fchk(_text("water_fk_chelpg.fchk"))
        self.assertIn("energy", reason)
        self.assertEqual(parser.data["mo_coeffs"], {})


class TestDamagedAndForeignFiles(unittest.TestCase):
    def test_incomplete_coordinate_array_is_rejected(self):
        sections = _module("fchk").read_sections(_text("water_fk_mo.fchk"))
        sections["Current cartesian coordinates"].pop()
        with self.assertRaisesRegex(ValueError, "coordinates"):
            _module("fchk").structure(sections)

    def test_other_text_is_not_a_checkpoint(self):
        fchk = _module("fchk")
        self.assertFalse(fchk.looks_like_fchk(_text("water_fk_mo.out")))
        self.assertFalse(fchk.looks_like_fchk(""))
        with self.assertRaises(ValueError):
            fchk.read_sections("just\nsome\ntext\n")

    def test_every_fixture_header_is_recognized_from_its_first_bytes(self):
        fchk = _module("fchk")
        self.assertTrue(fchk.looks_like_fchk(_text("water_fk_mo.fchk")[:512]))

    def test_a_truncated_array_is_reported(self):
        text = _text("water_fk_mo.fchk")
        cut = text.index("Alpha MO coefficients")
        with self.assertRaises(ValueError) as caught:
            _module("fchk").read_sections(text[: cut + 200])
        self.assertIn("Alpha MO coefficients", str(caught.exception))

    def test_text_sections_may_start_at_the_first_column(self):
        sections = _module("fchk").read_sections(_text("water_fk_chelpg.fchk"))
        self.assertEqual(sections["Route"], "#P HF/6-31G(d) Pop=ChelpG")
        self.assertEqual(sections["Gaussian Version"], "EM64W-G16RevB.01")

    def test_cartesian_g_is_marked_so_the_renderer_refuses_it(self):
        fchk = _module("fchk")
        shells = fchk.basis_shells(
            {
                "Shell types": [4],
                "Number of primitives per shell": [1],
                "Shell to atom map": [1],
                "Primitive exponents": [1.0],
                "Contraction coefficients": [1.0],
            },
            [[0.0, 0.0, 0.0]],
        )
        self.assertTrue(shells[0]["cartesian"])

    def test_the_esp_scheme_comes_from_the_route(self):
        fchk = _module("fchk")
        self.assertEqual(fchk._esp_label("#P HF/6-31G(d) Pop=ChelpG"), "CHELPG")
        self.assertEqual(fchk._esp_label("#P HF/6-31G(d) Pop=MK"), "Merz-Kollman")
        self.assertEqual(fchk._esp_label("#P HF/6-31G(d) Pop=(Full,MK)"), "Merz-Kollman")
        self.assertEqual(fchk._esp_label("#P HF/6-31G(d)"), "ESP")


if __name__ == "__main__":
    unittest.main()
