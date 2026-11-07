"""
tests/test_bond_analysis.py
Unit tests for the Bond Analysis panel's 3D-highlight logic and VDW scaling.

PySide6 and pyvista are stubbed so the module imports headlessly; the dialog's
Qt constructor is bypassed (object.__new__) to exercise the pure highlight
logic with a lightweight fake host (parser data + plotter).
"""

import os
import sys
import types
import importlib.util
import unittest
from unittest.mock import MagicMock


# ---------------------------------------------------------------------------
# Stub PySide6 + pyvista before importing bond_analysis
# ---------------------------------------------------------------------------


def _build_pv_stub():
    pv = types.ModuleType("pyvista")
    pv.Sphere = lambda *a, **k: object()

    class _LineObj:
        def tube(self, *a, **k):
            return object()

    pv.Line = lambda *a, **k: _LineObj()
    return pv


# A reusable pyvista stub; the highlight tests force it into sys.modules around
# each test because other test modules also (re)place "pyvista" in sys.modules.
_PV_STUB = _build_pv_stub()


def _install_stubs():
    for name in ["PySide6", "PySide6.QtWidgets", "PySide6.QtCore", "PySide6.QtGui"]:
        sys.modules.setdefault(name, types.ModuleType(name))

    class _Base:
        class DialogCode:
            Accepted = 1
            Rejected = 0

        def __init__(self, *a, **k):
            pass

        def __getattr__(self, n):
            return MagicMock()

    qtw = sys.modules["PySide6.QtWidgets"]
    for cls in [
        "QDialog",
        "QVBoxLayout",
        "QHBoxLayout",
        "QTabWidget",
        "QTableWidget",
        "QTableWidgetItem",
        "QHeaderView",
        "QLabel",
        "QMessageBox",
        "QPushButton",
        "QWidget",
    ]:
        setattr(qtw, cls, type(cls, (_Base,), {}))
    qtw.QApplication = MagicMock()
    sys.modules["PySide6.QtCore"].Qt = MagicMock()
    sys.modules["PySide6.QtGui"].QKeySequence = MagicMock()
    sys.modules["pyvista"] = _PV_STUB


_install_stubs()

_BOND_SRC = os.path.normpath(
    os.path.join(
        os.path.dirname(__file__), "..", "gaussian_result_analyzer", "bond_analysis.py"
    )
)
_spec = importlib.util.spec_from_file_location("bond_analysis_standalone", _BOND_SRC)
B = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(B)


class _Plotter:
    def __init__(self):
        self.added = []

    def add_mesh(self, *a, **k):
        obj = object()
        self.added.append(obj)
        return obj

    def remove_actor(self, actor):
        pass

    def render(self):
        pass


def _make_host(coords, atoms):
    return types.SimpleNamespace(
        parser=types.SimpleNamespace(data={"coords": coords, "atoms": atoms}),
        mw=types.SimpleNamespace(plotter=_Plotter()),
    )


def _bare_dialog(host):
    """A BondAnalysisDialog instance with the Qt constructor bypassed."""
    inst = B.BondAnalysisDialog.__new__(B.BondAnalysisDialog)
    inst._actors = []
    inst.parent_dlg = host
    return inst


# ---------------------------------------------------------------------------
# VDW radius scaling
# ---------------------------------------------------------------------------


class TestVdwScaling(unittest.TestCase):
    def test_fraction_is_40_percent(self):
        self.assertAlmostEqual(B._VDW_FRACTION, 0.40)

    def test_vdw_uses_periodic_table(self):
        saved = B._PERIODIC_TABLE
        try:
            B._PERIODIC_TABLE = MagicMock()
            B._PERIODIC_TABLE.GetRvdw.return_value = 1.52
            self.assertAlmostEqual(B._vdw("O"), 1.52)
            B._PERIODIC_TABLE.GetRvdw.assert_called_with("O")
        finally:
            B._PERIODIC_TABLE = saved

    def test_vdw_fallback_when_no_table(self):
        saved = B._PERIODIC_TABLE
        try:
            B._PERIODIC_TABLE = None
            self.assertAlmostEqual(B._vdw("O"), 1.70)
        finally:
            B._PERIODIC_TABLE = saved

    def test_halo_radius_is_fraction_times_vdw(self):
        saved = B._PERIODIC_TABLE
        try:
            B._PERIODIC_TABLE = None  # -> _vdw == 1.70
            inst = _bare_dialog(_make_host([[0, 0, 0]], ["O"]))
            self.assertAlmostEqual(inst._halo_radius(0, ["O"]), 0.40 * 1.70)
        finally:
            B._PERIODIC_TABLE = saved


# ---------------------------------------------------------------------------
# 3D highlight actor management
# ---------------------------------------------------------------------------


class TestHighlight(unittest.TestCase):
    def setUp(self):
        # Other test modules may have replaced sys.modules["pyvista"]; the
        # highlight methods lazily `import pyvista`, so pin our stub here.
        self._saved_pv = sys.modules.get("pyvista")
        sys.modules["pyvista"] = _PV_STUB
        self.host = _make_host([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], ["O", "H"])
        self.dlg = _bare_dialog(self.host)

    def tearDown(self):
        if self._saved_pv is not None:
            sys.modules["pyvista"] = self._saved_pv
        else:
            sys.modules.pop("pyvista", None)


    def test_atom_highlight_one_sphere_per_atom(self):
        self.dlg._highlight_atoms([0, 1])
        self.assertEqual(len(self.dlg._actors), 2)

    def test_atom_highlight_skips_out_of_range(self):
        self.dlg._highlight_atoms([0, 99])
        self.assertEqual(len(self.dlg._actors), 1)

    def test_bond_highlight_out_of_range_noop(self):
        self.dlg._highlight_bond(0, 99)
        self.assertEqual(len(self.dlg._actors), 0)


    def test_reselect_replaces_previous_highlight(self):
        self.dlg._highlight_bond(0, 1)
        self.dlg._highlight_atoms([0])  # clears bond first, then 1 sphere
        self.assertEqual(len(self.dlg._actors), 1)

    def test_no_plotter_is_safe(self):
        host = types.SimpleNamespace(
            parser=types.SimpleNamespace(data={"coords": [[0, 0, 0]], "atoms": ["O"]}),
            mw=types.SimpleNamespace(),  # no plotter attribute
        )
        dlg = _bare_dialog(host)
        dlg._highlight_atoms([0])  # must not raise
        self.assertEqual(len(dlg._actors), 0)


# ---------------------------------------------------------------------------
# Hybridization column text reflects the parsed (file) values
# ---------------------------------------------------------------------------


class TestHybridPercentColumn(unittest.TestCase):
    def test_bond_two_atoms_from_file_values(self):
        # O-H bond from nbo-test.out: O s22.56/p77.29/d0.15, H s99.76/p0.24
        hybrids = [
            {"atom_sym": "O", "s_pct": 22.56, "p_pct": 77.29, "d_pct": 0.15},
            {"atom_sym": "H", "s_pct": 99.76, "p_pct": 0.24, "d_pct": 0.0},
        ]
        # d < 1% is omitted; percentages are rounded to whole numbers.
        self.assertEqual(B._hyb_percent(hybrids), "O 23s 77p   ·   H 100s 0p")

    def test_d_shown_only_when_significant(self):
        hybrids = [{"atom_sym": "X", "s_pct": 20.0, "p_pct": 70.0, "d_pct": 10.0}]
        self.assertEqual(B._hyb_percent(hybrids), "X 20s 70p 10d")

    def test_empty_hybrids(self):
        self.assertEqual(B._hyb_percent([]), "")


if __name__ == "__main__":
    unittest.main()
