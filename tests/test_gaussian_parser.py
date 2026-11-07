"""GaussianParser against real Gaussian 16 output (tests/sample_outputs).

The expected numbers were read from the logs themselves. The MO test is the
independent one: coefficients and basis from the parser, evaluated by the engine
on a real-space grid, must give orthonormal orbitals, which fails if an order,
a normalization or a Cartesian/pure assignment is wrong.
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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gui_harness  # noqa: E402

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_SAMPLES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sample_outputs")
_PKG = "_gaussian_parser_only"
_BOHR = 0.529177249


def _load_module(name):
    """Import one package module without running the plugin's __init__."""
    if _PKG not in sys.modules:
        package = types.ModuleType(_PKG)
        package.__path__ = [os.path.join(_ROOT, "gaussian_result_analyzer")]
        package.PLUGIN_VERSION = "test"
        sys.modules[_PKG] = package
    return importlib.import_module(f"{_PKG}.{name}")


def _parse(name):
    parser = _load_module("parser").GaussianParser()
    with open(os.path.join(_SAMPLES, name + ".out"), encoding="utf-8") as handle:
        parser.load_from_memory(handle.read(), name + ".out")
    return parser.data


class TestBasics(unittest.TestCase):
    def test_version(self):
        self.assertEqual(_parse("water_mp2")["version"], "Gaussian 16 EM64W-G16RevB.01")

    def test_scf_energy(self):
        self.assertAlmostEqual(_parse("water_mp2")["scf_energy"], -76.0097717063)

    def test_geometry_is_in_the_standard_orientation(self):
        data = _parse("water_mp2")
        self.assertEqual(data["atoms"], ["O", "H", "H"])
        # Standard orientation of water: oxygen on the z axis.
        self.assertAlmostEqual(data["coords"][0][0], 0.0)
        self.assertAlmostEqual(data["coords"][0][1], 0.0)

    def test_charge_and_multiplicity(self):
        data = _parse("o2_triplet_uhf_opt")
        self.assertEqual((data["charge"], data["mult"]), (0, 3))

    def test_normal_termination(self):
        self.assertEqual(
            _parse("water_mp2")["termination_status"], "Terminated normally"
        )

    def test_error_termination(self):
        self.assertEqual(_parse("o2_triplet_opt")["termination_status"], "ERROR")


class TestTrajectories(unittest.TestCase):
    def test_optimization_cycles(self):
        data = _parse("benzene_opt_freq")
        self.assertEqual(len(data["scan_steps"]), 3)
        self.assertTrue(data["converged"])

    def test_relaxed_scan_keeps_only_the_converged_points(self):
        data = _parse("ethane_scan")
        self.assertTrue(data["is_scan"])
        self.assertEqual(len(data["scan_steps"]), 28)
        self.assertEqual({s["type"] for s in data["scan_steps"]}, {"scan_step"})

    def test_gradients_cover_every_atom(self):
        data = _parse("benzene_gradient")
        self.assertEqual(len(data["gradients"]), 12)
        # Without an external field the forces on a molecule sum to zero.
        total = [sum(g["vector"][k] for g in data["gradients"]) for k in range(3)]
        for value in total:
            self.assertAlmostEqual(value, 0.0, places=5)


class TestSpectra(unittest.TestCase):
    def test_frequency_count(self):
        self.assertEqual(len(_parse("benzene_opt_freq")["frequencies"]), 30)

    def test_raman_activity_is_read(self):
        modes = _parse("water_raman")["frequencies"]
        self.assertEqual(len(modes), 3)
        self.assertAlmostEqual(modes[0]["freq"], 1826.5605)
        self.assertAlmostEqual(modes[0]["ir"], 107.2699)
        self.assertAlmostEqual(modes[0]["raman"], 5.7237)

    def test_modes_have_one_vector_per_atom(self):
        for mode in _parse("benzene_opt_freq")["frequencies"]:
            self.assertEqual(len(mode["vector"]), 12)

    def test_thermochemistry_is_in_hartree(self):
        thermal = _parse("water_raman")["thermal"]
        self.assertAlmostEqual(thermal["temperature"], 298.15)
        self.assertAlmostEqual(thermal["zpe"], 0.022977, places=5)
        self.assertLess(thermal["gibbs"], thermal["enthalpy"])

    def test_tddft_states(self):
        states = _parse("h2co_tddft")["tddft"]
        self.assertEqual(len(states), 12)
        self.assertEqual(states[0]["label"], "Triplet-A2")
        self.assertAlmostEqual(states[0]["energy_ev"], 3.1664)

    def test_nmr_shieldings_and_couplings(self):
        data = _parse("methanol_nmr")
        self.assertEqual(len(data["nmr_shielding"]), 6)
        self.assertAlmostEqual(data["nmr_shielding"][0]["shielding"], 141.5334)
        self.assertEqual(len(data["nmr_couplings"]), 15)


class TestProperties(unittest.TestCase):
    def test_dipole(self):
        dipole = _parse("water_mp2")["dipole"]
        self.assertAlmostEqual(dipole["magnitude_debye"], 2.2443)

    def test_dispersion_energy(self):
        self.assertAlmostEqual(_parse("water_d3bj")["dispersion"], -0.0005738947)

    def test_spin_contamination(self):
        s2 = _parse("o2_triplet_uhf_opt")["spin_s2"]
        self.assertAlmostEqual(s2["actual"], 2.0066)
        self.assertAlmostEqual(s2["ideal"], 2.0)

    def test_correlation_energies(self):
        labels = {c["label"] for c in _parse("water_ccsdt")["energy_components"]}
        self.assertIn("MP2 Correlation", labels)
        self.assertEqual(len(labels), 4)

    def test_charge_schemes_sum_to_the_total_charge(self):
        expect = {
            "water_chelpg": "CHELPG",
            "water_mk": "Merz-Kollman",
            "water_raman": "APT",
            "water_6d_mo": "Hirshfeld",
        }
        for name, scheme in expect.items():
            with self.subTest(scheme=scheme):
                rows = _parse(name)["charges"][scheme]
                self.assertEqual(len(rows), 3)
                self.assertAlmostEqual(sum(r["charge"] for r in rows), 0.0, places=3)

    def test_nbo(self):
        data = _parse("water_6d_mo")
        self.assertEqual(len(data["nbo_orbitals"]), 19)
        self.assertEqual(len(data["nbo_perturbation"]), 4)


class TestOrbitals(unittest.TestCase):
    def test_restricted_mo_table_is_complete(self):
        data = _parse("benzene_opt_freq")
        self.assertEqual(len(data["mo_coeffs"]), 102)
        self.assertTrue(data["cartesian_d"])

    def test_unrestricted_table_is_complete_with_long_symmetry_labels(self):
        # (DLTG)-- has no room for the occupancy letter and must not end the table.
        data = _parse("o2_triplet_uhf_opt")
        self.assertEqual(len(data["mo_coeffs"]), 60)
        self.assertEqual(len(data["mo_coeffs"]["29_beta"]["coeffs"]), 30)

    def test_pure_shells_are_not_flagged_cartesian(self):
        data = _parse("h2co_5d7f_mo")
        self.assertFalse(data["cartesian_d"])
        self.assertFalse(data["cartesian_f"])


def _grid(centers, n_rad=90, n_ang=194):
    """Becke-partitioned atomic grid, in bohr."""
    theta = np.arccos(1 - 2 * (np.arange(n_ang) + 0.5) / n_ang)
    phi = np.pi * (1 + 5**0.5) * np.arange(n_ang)
    sphere = np.c_[
        np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)
    ]
    k = np.arange(1, n_rad + 1)
    x = np.cos(np.pi * k / (n_rad + 1))
    r = (1 + x) / (1 - x)
    wr = (
        (np.pi / (n_rad + 1))
        * np.sin(np.pi * k / (n_rad + 1)) ** 2
        * 2
        * r**2
        / (np.sqrt(1 - x**2) * (1 - x) ** 2)
    )
    points, weights = [], []
    for ia, center in enumerate(centers):
        p = center + (r[:, None, None] * sphere[None]).reshape(-1, 3)
        w = np.repeat(wr, n_ang) * 4 * np.pi / n_ang
        d = np.linalg.norm(p[:, None] - centers[None], axis=2)
        cell = np.ones((len(p), len(centers)))
        for i in range(len(centers)):
            for j in range(len(centers)):
                if i != j:
                    mu = (d[:, i] - d[:, j]) / np.linalg.norm(centers[i] - centers[j])
                    for _ in range(3):
                        mu = 1.5 * mu - 0.5 * mu**3
                    cell[:, i] *= 0.5 * (1 - mu)
        points.append(p)
        weights.append(w * cell[:, ia] / cell.sum(axis=1))
    return np.vstack(points), np.concatenate(weights)


@unittest.skipIf(np is None, "numpy is not installed")
class TestOrthonormalOnAGrid(unittest.TestCase):
    def _overlap(self, name, spin, count):
        engine_module = gui_harness.load_isolated("mo_engine")
        data = _parse(name)
        shells = [
            {
                "type": s["l"],
                "center": np.array(s["origin"]) / _BOHR,
                "exps": np.array(s["exps"]),
                "coeffs": np.array(s["coeffs"]),
                "cartesian": s["cartesian"],
            }
            for s in data["basis_set_shells"]
        ]
        engine = engine_module.BasisSetEngine(shells)
        points, weights = _grid(np.array(data["coords"]) / _BOHR)
        values = np.array(
            [
                engine.evaluate_mo_on_grid(
                    i,
                    points,
                    np.array(
                        [c["coeff"] for c in data["mo_coeffs"][f"{i}_{spin}"]["coeffs"]]
                    ),
                )
                for i in range(count)
            ]
        )
        return (values * weights) @ values.T

    def _assert_orthonormal(self, name, spin, count):
        overlap = self._overlap(name, spin, count)
        self.assertLess(np.abs(overlap - np.eye(count)).max(), 0.005)

    def test_cartesian_d(self):
        self._assert_orthonormal("water_6d_mo", "restricted", 19)

    def test_pure_d_and_f(self):
        self._assert_orthonormal("h2co_5d7f_mo", "restricted", 40)

    def test_unrestricted_alpha_and_beta(self):
        self._assert_orthonormal("o2_triplet_uhf_opt", "alpha", 30)
        self._assert_orthonormal("o2_triplet_uhf_opt", "beta", 30)


if __name__ == "__main__":
    unittest.main()
