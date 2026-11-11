"""Parser cases found by running a wider range of Gaussian 16 jobs through it.

Each fixture is a small job that shows the layout: a job that prints no MO table,
a stability analysis, a #T log, a double hybrid, a molecule above 50 atoms,
a rigid scan, an IRC, CASSCF, counterpoise and anharmonic jobs.
"""

import importlib
import os
import sys
import types
import unittest

try:
    import numpy as np
except ImportError:
    np = None

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_SAMPLES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sample_outputs")
_PKG = "_gaussian_parser_cases"


def _module(name):
    if _PKG not in sys.modules:
        package = types.ModuleType(_PKG)
        package.__path__ = [os.path.join(_ROOT, "gaussian_result_analyzer")]
        package.PLUGIN_VERSION = "test"
        sys.modules[_PKG] = package
    return importlib.import_module(f"{_PKG}.{name}")


def _parse(name):
    parser = _module("parser").GaussianParser()
    with open(os.path.join(_SAMPLES, name + ".out"), encoding="utf-8") as handle:
        parser.load_from_memory(handle.read(), name + ".out")
    return parser.data


class TestEigenvalueSummary(unittest.TestCase):
    """Without Pop=Full the energies come from the occ./virt. eigenvalue lines."""

    def test_every_orbital_is_read(self):
        data = _parse("water_polar")
        self.assertEqual(data["mo_coeffs"], {})
        energies = [o["energy"] for o in data["orbital_energies"]]
        self.assertEqual(len(energies), 19)
        self.assertEqual(energies, sorted(energies))
        self.assertAlmostEqual(energies[0], -19.13605)
        self.assertAlmostEqual(energies[4], -0.29182)
        self.assertAlmostEqual(energies[5], 0.06564)
        self.assertAlmostEqual(energies[18], 3.55828)

    def test_occupations_follow_the_electron_count(self):
        orbitals = _parse("water_polar")["orbital_energies"]
        self.assertEqual([o["occ"] for o in orbitals[:6]], [2.0] * 5 + [0.0])

    def test_they_match_the_checkpoint_of_a_like_job(self):
        log = _parse("water_fk_plain")["orbital_energies"]
        parser = _module("parser").GaussianParser()
        with open(os.path.join(_SAMPLES, "water_fk_plain.fchk"), encoding="utf-8") as f:
            parser.load_fchk(f.read(), "water_fk_plain.fchk")
        fchk = parser.data["orbital_energies"]
        self.assertEqual(len(log), len(fchk))
        for a, b in zip(log, fchk):
            self.assertAlmostEqual(a["energy"], b["energy"], delta=6e-6)
            self.assertEqual(a["occ"], b["occ"])

    def test_unrestricted_summary(self):
        orbitals = _parse("methyl_radical_nbo")["orbital_energies"]
        alpha = [o for o in orbitals if o["spin"] == "alpha"]
        beta = [o for o in orbitals if o["spin"] == "beta"]
        self.assertEqual((len(alpha), len(beta)), (21, 21))
        self.assertEqual(sum(o["occ"] for o in alpha), 5.0)
        self.assertEqual(sum(o["occ"] for o in beta), 4.0)


class TestSpinContamination(unittest.TestCase):
    def test_a_stability_analysis_of_a_singlet_reports_none(self):
        """Its eigenvectors carry <S**2>=2.000, which is not the wavefunction's."""
        self.assertIsNone(_parse("water_stable")["spin_s2"])

    def test_a_radical_reports_the_converged_value(self):
        s2 = _parse("methyl_radical_nbo")["spin_s2"]
        self.assertAlmostEqual(s2["actual"], 0.7538, places=3)
        self.assertAlmostEqual(s2["ideal"], 0.75)


class TestTerseLog(unittest.TestCase):
    """#T prints the thermochemistry without its title line."""

    def test_the_fixture_has_no_title(self):
        with open(os.path.join(_SAMPLES, "water_terse_t.out"), encoding="utf-8") as f:
            text = f.read()
        self.assertNotIn("- Thermochemistry -", text)
        self.assertIn("Zero-point correction=", text)

    def test_thermochemistry_is_read(self):
        thermal = _parse("water_terse_t")["thermal"]
        self.assertAlmostEqual(thermal["zpe"], 0.021163)
        self.assertAlmostEqual(thermal["gibbs"], -76.405458)
        self.assertAlmostEqual(thermal["temperature"], 298.15)

    def test_frequencies_are_read(self):
        self.assertEqual(len(_parse("water_terse_t")["frequencies"]), 3)


class TestDoubleHybrid(unittest.TestCase):
    """The total is on the E2(B2PLYP) line; SCF Done lacks the PT2 part."""

    def test_total_energy(self):
        rows = {r["label"]: r["value"] for r in _parse("water_b2plyp")["energy_components"]}
        self.assertAlmostEqual(rows["Double-Hybrid Total"], -76.334890892898)
        self.assertAlmostEqual(rows["Double-Hybrid Correlation"], -0.05997931296)

    def test_other_jobs_have_no_such_row(self):
        labels = [r["label"] for r in _parse("water_mp2")["energy_components"]]
        self.assertNotIn("Double-Hybrid Total", labels)


#: Rows of the job's formatted checkpoint, which stores the gradient in the
#: standard orientation.
EXPECTED_GRADIENT = {0: [0.02431704, -9.739e-05, 0.00226528], 16: [0.01510906, -0.01795279, 0.01063452], 17: [-0.00271396, -0.0052722, 0.0043508], 52: [0.00214419, 0.00014741, -0.00765234]}


@unittest.skipIf(np is None, "numpy is required to rotate the forces")
class TestForcesAboveFiftyAtoms(unittest.TestCase):
    """Gaussian then prints no input orientation; the input section gives the frame."""

    def setUp(self):
        self.data = _parse("heptadecane_force")

    def test_the_fixture_prints_only_the_standard_orientation(self):
        with open(os.path.join(_SAMPLES, "heptadecane_force.out"), encoding="utf-8") as f:
            self.assertNotIn("Input orientation:", f.read())
        self.assertEqual(len(self.data["atoms"]), 53)

    def test_the_gradient_is_in_the_frame_of_the_structure(self):
        for index, expected in EXPECTED_GRADIENT.items():
            for got, want in zip(self.data["gradients"][index]["vector"], expected):
                self.assertAlmostEqual(got, want, places=6)

    def test_net_torque_vanishes_in_that_frame(self):
        coords = np.array(self.data["coords"])
        grads = np.array([g["vector"] for g in self.data["gradients"]])
        torque = np.cross(coords - coords.mean(axis=0), grads).sum(axis=0)
        self.assertLess(np.abs(torque).max(), 1e-4)


#: Rows of the 104-atom optimization's checkpoint gradient (the last step).
EXPECTED_OPT_GRADIENT = {
    0: [0.00128125, -0.00080145, -0.00055262],
    33: [-0.00103036, 0.00107966, 0.00089364],
    103: [0.00016096, 0.00036979, -0.00055571],
}


@unittest.skipIf(np is None, "numpy is required to rotate the forces")
class TestOptimizationAboveFiftyAtoms(unittest.TestCase):
    """Every step, not only a single-point one, is turned into the structure's frame."""

    def setUp(self):
        self.data = _parse("tetratriacontane_opt")

    def test_the_final_gradient_matches_the_checkpoint(self):
        self.assertEqual(len(self.data["atoms"]), 104)
        for index, expected in EXPECTED_OPT_GRADIENT.items():
            for got, want in zip(self.data["gradients"][index]["vector"], expected):
                self.assertAlmostEqual(got, want, places=7)

    def test_every_step_has_torque_free_forces(self):
        steps = self.data["scan_steps"]
        self.assertEqual(len(steps), 3)
        for step in steps:
            coords = np.array(step["coords"])
            grads = np.array([g["vector"] for g in step["gradients"]])
            torque = np.cross(coords - coords.mean(axis=0), grads).sum(axis=0)
            self.assertLess(np.abs(torque).max(), 1e-6)


class TestRigidScan(unittest.TestCase):
    def test_one_point_per_energy_with_its_value(self):
        data = _parse("ethane_rigid_scan")
        self.assertTrue(data["is_scan"])
        points = [(s["scan_coord"], s["energy"]) for s in data["scan_steps"]]
        self.assertEqual(
            points,
            [
                (60.0, -78.7935626917),
                (80.0, -78.792307453),
                (100.0, -78.7897872295),
                (120.0, -78.7885222839),
            ],
        )
        self.assertTrue(all(len(s["atoms"]) == 8 for s in data["scan_steps"]))


class TestIRC(unittest.TestCase):
    """The checkpoint stores 13 points: 6 each way and the TS."""

    def setUp(self):
        self.data = _parse("hcn_irc")
        self.steps = self.data["scan_steps"]

    def test_the_path_runs_from_the_reverse_to_the_forward_end(self):
        coords = [s["scan_coord"] for s in self.steps]
        self.assertEqual(len(coords), 13)
        self.assertEqual(coords, sorted(coords))
        self.assertAlmostEqual(coords[0], -0.64901)
        self.assertEqual(coords[6], 0.0)
        self.assertAlmostEqual(coords[-1], 0.64896)

    def test_the_ts_is_the_highest_point(self):
        energies = [s["energy"] for s in self.steps]
        self.assertAlmostEqual(energies[6], -92.791953409)
        self.assertEqual(max(energies), energies[6])
        self.assertAlmostEqual(energies[0], -92.8017964484)

    def test_the_final_structure_is_the_last_point_computed(self):
        self.assertEqual(self.data["coords"][2], [0.079756, 1.03828, 0.0])

    def test_the_ts_search_before_it_is_not_part_of_the_path(self):
        self.assertTrue(all(s["type"] == "scan_step" for s in self.steps))


class TestJobsWithAuxiliaryCalculations(unittest.TestCase):
    def test_counterpoise_reports_the_full_system(self):
        data = _parse("water_dimer_cp")
        self.assertAlmostEqual(data["scf_energy"], -152.814919518)
        self.assertEqual(len(data["orbital_energies"]), 38)
        rows = {r["label"]: r["value"] for r in data["energy_components"]}
        self.assertAlmostEqual(rows["Counterpoise Corrected"], -152.813174427068)
        self.assertAlmostEqual(rows["BSSE"], 0.001745090539)

    def test_anharmonic_reports_the_equilibrium_energy(self):
        data = _parse("water_freq_anharm")
        self.assertAlmostEqual(data["scf_energy"], -76.0107464985)
        self.assertEqual(len(data["frequencies"]), 3)

    def test_casscf_energy(self):
        self.assertAlmostEqual(_parse("n2_casscf")["scf_energy"], -109.0787529468)


if __name__ == "__main__":
    unittest.main()
