"""
tests/test_traj_external.py
Covers TrajectoryResultDialog.load_external_trj: replacing the in-memory
trajectory with frames read from an XYZ file, and the view state that has to
be re-derived afterwards (scan-point toggles, coordinate axis availability,
slider extent and the energy baseline).
"""

import os
import sys
import tempfile
import unittest
from unittest.mock import MagicMock, patch

sys.path.insert(0, os.path.dirname(__file__))
import gui_harness  # noqa: E402

TR = gui_harness.load_isolated("traj_analysis")


def _xyz(frames):
    """Multi-frame XYZ; the comment line carries the energy, as Gaussian writes it."""
    out = []
    for energy, x in frames:
        out.append("2")
        out.append(f"Coordinates from Gaussian-job E {energy:.8f}")
        out.append("O    0.000000    0.000000    0.000000")
        out.append(f"H    {x:.6f}    0.000000    0.000000")
    return "\n".join(out) + "\n"


def _step(energy, x=0.0, stype="opt_cycle", **extra):
    s = {
        "energy": energy,
        "atoms": ["O", "H"],
        "coords": [[0.0, 0.0, 0.0], [x, 0.0, 0.0]],
        "type": stype,
    }
    s.update(extra)
    return s


class _ExternalCase(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name
        self.addCleanup(self._tmp.cleanup)

        self.gl = MagicMock()
        self.dlg = TR.TrajectoryResultDialog(
            self.gl,
            [_step(-100.0, 0.96), _step(-100.5, 0.97)],
            context=MagicMock(),
        )

    def _write(self, frames, name="traj.xyz"):
        path = os.path.join(self.tmp, name)
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(_xyz(frames))
        return path


class TestLoadExternalTrajectory(_ExternalCase):
    def test_frames_replace_the_in_memory_steps(self):
        path = self._write([(-100.0, 0.96), (-100.4, 0.97), (-100.6, 0.98)])
        self.dlg.load_external_trj(path)
        self.assertEqual(len(self.dlg.steps), 3)

    def test_geometry_comes_from_the_file(self):
        path = self._write([(-100.0, 0.96), (-100.4, 1.20)])
        self.dlg.load_external_trj(path)
        self.assertAlmostEqual(self.dlg.steps[1]["coords"][1][0], 1.20, places=5)

    def test_atoms_come_from_the_file(self):
        path = self._write([(-100.0, 0.96)])
        self.dlg.load_external_trj(path)
        self.assertEqual(self.dlg.steps[0]["atoms"], ["O", "H"])

    def test_zero_energy_frames_are_discarded(self):
        path = self._write([(-100.0, 0.96), (0.0, 0.97), (-100.6, 0.98)])
        self.dlg.load_external_trj(path)
        self.assertEqual(len(self.dlg.steps), 2)

    def test_the_energy_baseline_is_recomputed(self):
        path = self._write([(-200.0, 0.96), (-201.0, 0.97)])
        self.dlg.load_external_trj(path)
        self.assertAlmostEqual(self.dlg.global_min_e, -201.0)

    def test_the_full_history_is_kept_alongside_the_scan_points(self):
        path = self._write([(-100.0, 0.96), (-100.4, 0.97)])
        self.dlg.load_external_trj(path)
        self.assertEqual(len(self.dlg.all_steps), 2)

    def test_without_scan_points_the_view_toggle_is_disabled(self):
        path = self._write([(-100.0, 0.96), (-100.4, 0.97)])
        self.dlg.load_external_trj(path)
        # every frame is a scan point, so there is nothing to toggle between
        self.assertFalse(self.dlg.radio_scan.isEnabled())
        self.assertFalse(self.dlg.showing_scan_points)

    def test_the_slider_spans_the_new_frames(self):
        path = self._write([(-100.0, 0.96), (-100.4, 0.97), (-100.6, 0.98)])
        self.dlg.slider = MagicMock()
        self.dlg.load_external_trj(path)
        self.dlg.slider.setMaximum.assert_called_with(2)

    def test_playback_is_enabled_once_structures_are_present(self):
        path = self._write([(-100.0, 0.96), (-100.4, 0.97)])
        self.dlg.load_external_trj(path)
        self.assertTrue(self.dlg.btn_play.isEnabled())

    def test_external_frames_always_offer_the_coordinate_axis(self):
        # parse_xyz_content tags every frame "neb_step", and a path is
        # coordinate-based even when the XYZ comments carry no distance.
        path = self._write([(-100.0, 0.96), (-100.4, 0.97)])
        self.dlg.load_external_trj(path)
        self.assertTrue(self.dlg.show_coord_x)
        self.assertTrue(self.dlg.radio_coord.isEnabled())

    def test_path_distances_are_carried_over_when_the_frame_count_matches(self):
        # An NEB path summary knows the reaction coordinate; the XYZ does not.
        self.dlg.steps = [
            _step(-100.0, 0.96, "neb_image", dist=0.0),
            _step(-100.5, 0.97, "neb_image", dist=0.5),
        ]
        path = self._write([(-100.0, 0.96), (-100.5, 0.97)])
        self.dlg.load_external_trj(path)
        self.assertEqual([s.get("dist") for s in self.dlg.steps], [0.0, 0.5])

    def test_carried_over_distances_select_the_coordinate_axis(self):
        self.dlg.steps = [
            _step(-100.0, 0.96, "neb_image", dist=0.0),
            _step(-100.5, 0.97, "neb_image", dist=0.5),
        ]
        path = self._write([(-100.0, 0.96), (-100.5, 0.97)])
        self.dlg.load_external_trj(path)
        self.assertTrue(self.dlg.show_coord_x)

    def test_distances_are_not_carried_over_when_the_counts_differ(self):
        self.dlg.steps = [_step(-100.0, 0.96, "neb_image", dist=0.0)]
        path = self._write([(-100.0, 0.96), (-100.5, 0.97)])
        self.dlg.load_external_trj(path)
        self.assertEqual([s.get("dist") for s in self.dlg.steps], [None, None])

    def test_the_viewer_is_told_the_geometry_is_fixed(self):
        path = self._write([(-100.0, 0.96)])
        self.dlg.load_external_trj(path)
        self.assertTrue(self.gl.is_xyz_derived)

    def test_the_viewer_is_switched_into_3d_mode(self):
        path = self._write([(-100.0, 0.96)])
        self.dlg.load_external_trj(path)
        self.gl.ui_manager.enter_3d_viewer_mode.assert_called_once()

    def test_an_empty_file_is_reported(self):
        path = os.path.join(self.tmp, "empty.xyz")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("")
        before = list(self.dlg.steps)
        with patch.object(TR.QMessageBox, "warning") as warn:
            self.dlg.load_external_trj(path)
        warn.assert_called_once()
        self.assertEqual(self.dlg.steps, before)

    def test_an_empty_file_can_be_rejected_quietly(self):
        path = os.path.join(self.tmp, "empty.xyz")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("")
        with patch.object(TR.QMessageBox, "warning") as warn:
            self.dlg.load_external_trj(path, silent=True)
        warn.assert_not_called()

    def test_a_missing_file_leaves_the_trajectory_alone(self):
        before = list(self.dlg.steps)
        self.dlg.load_external_trj(os.path.join(self.tmp, "nope.xyz"), silent=True)
        self.assertEqual(self.dlg.steps, before)


if __name__ == "__main__":
    unittest.main()
