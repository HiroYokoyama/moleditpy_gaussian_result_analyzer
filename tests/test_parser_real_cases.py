"""Parser cases found by running a wider range of Gaussian 16 jobs through it.

Each fixture is a small job that shows the layout: a job that prints no MO table,
a stability analysis, a #T log, a double hybrid, a molecule above 50 atoms.
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


if __name__ == "__main__":
    unittest.main()
