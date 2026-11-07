"""
tests/test_parser_xyz_regressions.py
Regression tests for parse_xyz_content bugs fixed in 3.8.5, plus extra
coverage of the comment-line energy/coordinate extraction.

Bug 1: the energy-label regex had no word boundary, so the trailing "e" of a
preceding word (e.g. "Coordinate 1.2 Energy -3.4") matched the bare-"E"
alternative case-insensitively and stole the wrong number as the energy.

Bug 2: a single malformed coordinate line raised ValueError and aborted the
entire multi-frame parse instead of skipping that atom.

parser.py imports only stdlib — GaussianParser is loaded directly without stubs.
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(__file__))
from _parser_loader import load_standalone_parser  # noqa: E402


def _load_parser():
    return load_standalone_parser("gaussian_parser_xyz_regression")


GaussianParser = _load_parser().GaussianParser


def _frame(comment, atom_lines=("O 0.0 0.0 0.11779", "H 0.0 0.75545 -0.47116")):
    return "\n".join([str(len(atom_lines)), comment, *atom_lines]) + "\n"


class TestXYZFrameValidation(unittest.TestCase):
    def test_invalid_counts_cannot_rewind_the_parser(self):
        class BoundedLines(list):
            reads = 0

            def __getitem__(self, index):
                self.reads += 1
                if self.reads > 100:
                    raise AssertionError("XYZ parser did not advance")
                return super().__getitem__(index)

        class BoundedText(str):
            def splitlines(self, keepends=False):
                return BoundedLines(super().splitlines(keepends))

        content = BoundedText("-2\nTS\n-1\nCI\n0\nCI\n" + _frame("E -76.1"))
        steps = GaussianParser().parse_xyz_content(content)
        self.assertEqual(len(steps), 1)
        self.assertEqual(steps[0]["atoms"], ["O", "H"])

    def test_scientific_notation_in_labels(self):
        steps = GaussianParser().parse_xyz_content(_frame("Energy: -1.23e+02 Dist 1e-2"))
        self.assertAlmostEqual(steps[0]["energy"], -123.0)
        self.assertAlmostEqual(steps[0]["scan_coord"], 0.01)

    def test_scientific_notation_in_fallback_energy(self):
        steps = GaussianParser().parse_xyz_content(_frame("Frame 1 -1.23E+02"))
        self.assertAlmostEqual(steps[0]["energy"], -123.0)


class TestEnergyLabelWordBoundary(unittest.TestCase):
    """Bug 1: bare-'E' alternative must not match the tail of another word."""

    def test_coordinate_before_energy_label(self):
        # Old behavior: 'coordinatE 1.2' matched first -> energy = 1.2
        steps = GaussianParser().parse_xyz_content(_frame("Coordinate 1.2 Energy -76.400"))
        self.assertAlmostEqual(steps[0]["energy"], -76.400, places=3)

    def test_word_ending_in_e_before_number_falls_back_to_last_float(self):
        # 'Frame 3 -76.5': no real energy label -> fallback = last float
        steps = GaussianParser().parse_xyz_content(_frame("Frame 3 -76.5"))
        self.assertAlmostEqual(steps[0]["energy"], -76.5, places=3)

    def test_bare_e_label_still_works(self):
        steps = GaussianParser().parse_xyz_content(_frame("Iter 5 E -76.345"))
        self.assertAlmostEqual(steps[0]["energy"], -76.345, places=3)

    def test_e_colon_label_still_works(self):
        steps = GaussianParser().parse_xyz_content(_frame("E: -76.123"))
        self.assertAlmostEqual(steps[0]["energy"], -76.123, places=3)

    def test_e_equals_label_still_works(self):
        steps = GaussianParser().parse_xyz_content(_frame("E=-76.222"))
        self.assertAlmostEqual(steps[0]["energy"], -76.222, places=3)

    def test_energy_label_case_insensitive(self):
        steps = GaussianParser().parse_xyz_content(_frame("energy: -76.888"))
        self.assertAlmostEqual(steps[0]["energy"], -76.888, places=3)


class TestDistLabelExtraction(unittest.TestCase):
    def test_dist_label(self):
        steps = GaussianParser().parse_xyz_content(_frame("Dist 1.25 E -76.1"))
        self.assertAlmostEqual(steps[0]["dist"], 1.25, places=3)
        self.assertAlmostEqual(steps[0]["scan_coord"], 1.25, places=3)

    def test_distance_label(self):
        steps = GaussianParser().parse_xyz_content(_frame("Distance: 2.50 E -76.1"))
        self.assertAlmostEqual(steps[0]["dist"], 2.50, places=3)

    def test_coord_label(self):
        steps = GaussianParser().parse_xyz_content(_frame("Coord 0.95 E -76.1"))
        self.assertAlmostEqual(steps[0]["dist"], 0.95, places=3)

    def test_scan_label(self):
        steps = GaussianParser().parse_xyz_content(_frame("Scan 1.80 E -76.1"))
        self.assertAlmostEqual(steps[0]["dist"], 1.80, places=3)

    def test_no_dist_label_is_none(self):
        steps = GaussianParser().parse_xyz_content(_frame("E -76.1"))
        self.assertIsNone(steps[0]["dist"])

    def test_word_containing_dist_prefix_not_matched(self):
        # Leading \b: 'undistorted' must not feed the Dist label.
        steps = GaussianParser().parse_xyz_content(_frame("undistorted E -76.1"))
        self.assertIsNone(steps[0]["dist"])


class TestMalformedCoordinateLine(unittest.TestCase):
    """Bug 2: one bad coordinate line must not abort the whole parse."""

    def test_bad_coord_line_skipped(self):
        content = _frame(
            "E -76.1",
            atom_lines=(
                "O 0.0 0.0 0.11779",
                "H 0.0 abc -0.47116",  # malformed y
                "H 0.0 -0.75545 -0.47116",
            ),
        )
        steps = GaussianParser().parse_xyz_content(content)
        self.assertEqual(len(steps), 1)
        self.assertEqual(steps[0]["atoms"], ["O", "H"])
        self.assertEqual(len(steps[0]["coords"]), 2)

    def test_bad_frame_does_not_kill_later_frames(self):
        content = _frame("E -76.1", atom_lines=("O 0.0 0.0 0.1", "H x y z")) + _frame(
            "E -76.2"
        )
        steps = GaussianParser().parse_xyz_content(content)
        self.assertEqual(len(steps), 2)
        self.assertAlmostEqual(steps[1]["energy"], -76.2, places=3)
        self.assertEqual(len(steps[1]["coords"]), 2)

    def test_atoms_and_coords_stay_aligned(self):
        content = _frame(
            "E -76.1",
            atom_lines=("O 0.0 0.0 0.1", "H bad bad bad", "H 1.0 2.0 3.0"),
        )
        steps = GaussianParser().parse_xyz_content(content)
        self.assertEqual(len(steps[0]["atoms"]), len(steps[0]["coords"]))
        self.assertAlmostEqual(steps[0]["coords"][1][2], 3.0, places=3)


if __name__ == "__main__":
    unittest.main()
