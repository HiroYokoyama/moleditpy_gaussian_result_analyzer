"""Layouts the fixtures do not show, made by editing a real log in memory."""

import importlib
import os
import sys
import types
import unittest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_SAMPLES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sample_outputs")
_PKG = "_gaussian_parser_audit"

_HPMODES = """ Frequencies ---  1826.6177  4070.4203  4188.6017
 Reduced masses ---     1.0823     1.0455     1.0828
 Force constants ---     2.1276    10.2058    11.1929
 IR Intensities ---   107.2795    18.1946    58.0142
 Coord Atom Element:
   1     1     8          0.00000   0.00000   0.00000
   2     1     8          0.00000   0.00000   0.00000
   3     1     8          0.06841   0.04942   0.00000
"""


def _module(name):
    if _PKG not in sys.modules:
        package = types.ModuleType(_PKG)
        package.__path__ = [os.path.join(_ROOT, "gaussian_result_analyzer")]
        package.PLUGIN_VERSION = "test"
        sys.modules[_PKG] = package
    return importlib.import_module(f"{_PKG}.{name}")


def _text(name):
    with open(os.path.join(_SAMPLES, name + ".out"), encoding="utf-8") as handle:
        return handle.read()


def _parse_text(text):
    parser = _module("parser").GaussianParser()
    parser.load_from_memory(text, "edited.out")
    return parser.data


class TestHPModes(unittest.TestCase):
    """Freq=HPModes prints "Frequencies ---" rows above the usual table."""

    def test_the_high_precision_table_does_not_stop_the_load(self):
        text = _text("water_raman")
        section = text.rfind("Harmonic frequencies (cm**-1)")
        first = text.find(" Frequencies --", section)
        plain = _parse_text(text)["frequencies"]
        data = _parse_text(text[:first] + _HPMODES + text[first:])
        self.assertEqual(
            [m["freq"] for m in data["frequencies"]], [m["freq"] for m in plain]
        )
        self.assertEqual(
            [m["vector"] for m in data["frequencies"]], [m["vector"] for m in plain]
        )


class TestChargeOfTheWholeSystem(unittest.TestCase):
    def test_a_counterpoise_fragment_charge_is_not_the_job_charge(self):
        text = _text("water_dimer_cp").replace(
            " Charge =  0 Multiplicity = 1 in fragment      2.",
            " Charge = -1 Multiplicity = 2 in fragment      2.",
        )
        data = _parse_text(text)
        self.assertEqual((data["charge"], data["mult"]), (0, 1))


class TestOlderChargeTitles(unittest.TestCase):
    def test_gaussian_09_titles(self):
        text = _text("water_raman")
        new = _parse_text(text)["charges"]
        old = _parse_text(
            text.replace(" Mulliken charges:", " Mulliken atomic charges:").replace(
                " APT charges:", " APT atomic charges:"
            )
        )["charges"]
        self.assertIn("Mulliken", new)
        self.assertIn("APT", new)
        self.assertEqual(old, new)


class TestWrappedRoute(unittest.TestCase):
    def test_a_keyword_broken_across_rows_is_read(self):
        text = _text("water_chelpg").replace(
            " #P HF/6-31G(d) Pop=ChelpG\n", " #P HF/6-31G(d) Pop=Che\n lpG\n", 1
        )
        self.assertIn("CHELPG", _parse_text(text)["charges"])

    def test_fopt_is_an_optimization(self):
        text = _text("water_fk_opt")
        steps = len(_parse_text(text)["scan_steps"])
        self.assertGreater(steps, 1)
        route = next(l for l in text.splitlines() if l.startswith(" #"))
        self.assertRegex(route, r"(?i)\bopt\b")
        edited = text.replace(route, route.replace("Opt", "FOpt").replace("opt", "fopt"), 1)
        self.assertEqual(len(_parse_text(edited)["scan_steps"]), steps)


class TestUnfinishedLaterStep(unittest.TestCase):
    """The first step's Normal termination says nothing about the second."""

    def test_an_internal_step_that_has_not_ended(self):
        lines = _text("water_raman").splitlines()
        step = next(i for i, l in enumerate(lines) if "Proceeding to internal job step" in l)
        data = _parse_text("\n".join(lines[: step + 40]))
        self.assertEqual(data["termination_status"], "Running")

    def test_a_link1_job_that_has_not_ended(self):
        lines = _text("hcn_irc").splitlines()
        ends = [i for i, l in enumerate(lines) if "Normal termination" in l]
        data = _parse_text("\n".join(lines[: ends[1] + 120]))
        self.assertEqual(data["termination_status"], "Running")

    def test_finished_jobs_are_unchanged(self):
        for name in ("water_raman", "hcn_irc"):
            self.assertEqual(
                _parse_text(_text(name))["termination_status"], "Terminated normally"
            )


if __name__ == "__main__":
    unittest.main()
