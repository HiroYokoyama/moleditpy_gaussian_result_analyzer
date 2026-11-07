"""Reading and parsing Gaussian output behind a progress dialog.

A large ``.out`` is read and parsed on the GUI thread, so without feedback the
application looks hung.
"""

import logging
import os

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication, QProgressDialog

from .fchk import FCHK_EXTENSIONS, looks_like_fchk
from .parser import GaussianParser, ParseCancelled

#: Tried in order for BOM-less files. UTF-16 is deliberately absent: almost any
#: even-length byte string decodes as UTF-16, so a cp1252 file used to come back
#: as CJK garbage instead of reaching the 8-bit codecs. latin-1 never fails, so
#: it must stay last.
ENCODINGS = ("utf-8", "cp1252", "latin-1")

_UTF16_BOMS = (b"\xff\xfe", b"\xfe\xff")

#: Reading counts as this fraction of the whole load in the progress bar.
_READ_FRACTION = 0.25

#: Don't flash a dialog for a small file that loads instantly (ms).
_MIN_DURATION_MS = 300


def read_gaussian_text(path: str) -> str:
    """Return the text of *path*, trying several encodings.

    UTF-16 is used only when the file starts with its byte-order mark.
    Undecodable bytes are replaced rather than failing the load; I/O errors
    propagate.
    """
    with open(path, "rb") as f:
        raw = f.read()
    if raw.startswith(_UTF16_BOMS):
        return raw.decode("utf-16", errors="replace").replace("\r\n", "\n")
    for enc in ENCODINGS:
        try:
            text = raw.decode(enc)
        except UnicodeError:
            continue
        return text.replace("\r\n", "\n").replace("\r", "\n")
    return raw.decode("utf-8", errors="replace")


class LoadProgress:
    """Modal progress dialog for one file load."""

    def __init__(self, parent, path):
        self._cancelled = False
        self._dlg = QProgressDialog(
            f"Reading {os.path.basename(path)}...", "Cancel", 0, 100, parent
        )
        self._dlg.setWindowTitle("Loading Gaussian Output")
        self._dlg.setWindowModality(Qt.WindowModality.WindowModal)
        self._dlg.setMinimumDuration(_MIN_DURATION_MS)
        self._dlg.setAutoClose(False)
        self._dlg.setAutoReset(False)
        # A signal, not wasCanceled(): polling a permissive test stub would
        # read as cancelled and abort every load.
        try:
            self._dlg.canceled.connect(self._on_cancel)
        except (AttributeError, TypeError) as exc:
            logging.debug("LoadProgress: cancel not wired — %s", exc)
        self._dlg.setValue(0)
        QApplication.processEvents()

    @property
    def cancelled(self):
        """Whether the user cancelled the progress dialog."""
        return self._cancelled

    def _on_cancel(self):
        """Record that the dialog's Cancel button was clicked."""
        self._cancelled = True

    def _pump(self, value, label):
        """Update the dialog's label and progress value and flush Qt events."""
        self._dlg.setLabelText(label)
        self._dlg.setValue(int(value))
        QApplication.processEvents()

    def read_done(self, path):
        """Advance the bar to the read/parse boundary once the file has been read."""
        self._pump(_READ_FRACTION * 100, f"Parsing {os.path.basename(path)}...")

    def on_parse_step(self, done, total, label):
        """Advance the bar for one parse_all step; raise if the user cancelled."""
        if self._cancelled:
            raise ParseCancelled(label)
        frac = _READ_FRACTION + (1.0 - _READ_FRACTION) * (done / max(total, 1))
        self._pump(frac * 100, label if done >= total else f"{label}...")

    def close(self):
        """Hide and dispose of the progress dialog."""
        try:
            self._dlg.hide()
            self._dlg.close()
            self._dlg.deleteLater()
        except (AttributeError, RuntimeError) as exc:
            logging.debug("LoadProgress: hide/close/delete failed — %s", exc)


def is_fchk_path(path):
    """Whether *path* has an extension formchk writes."""
    return path.lower().endswith(FCHK_EXTENSIONS)


def _first_existing(path, extensions):
    """The file beside *path* with one of *extensions*, in either case."""
    stem = os.path.splitext(path)[0]
    for extension in extensions:
        for candidate in (stem + extension, stem + extension.upper()):
            if os.path.isfile(candidate):
                return candidate
    return ""


def sibling_fchk(path):
    """The checkpoint with the same name as an output file, or an empty string."""
    return _first_existing(path, FCHK_EXTENSIONS)


def attach_fchk(parser, fchk_path):
    """Use a checkpoint's orbitals, recording why not when it is another job's.

    True when they were taken. A broken file beside the output must not stop the
    output from opening, so that outcome is a note rather than an exception.
    """
    name = os.path.basename(fchk_path)
    try:
        reason = parser.merge_fchk(read_gaussian_text(fchk_path))
    except (OSError, ValueError) as exc:
        logging.warning("Reading the checkpoint %s failed: %s", fchk_path, exc)
        parser.data["fchk_note"] = f"{name} could not be read: {exc}"
        return False
    if reason:
        logging.warning("Not using the checkpoint %s: %s", fchk_path, reason)
        parser.data["fchk_note"] = f"{name} not used: {reason}."
        return False
    parser.data["fchk_path"] = fchk_path
    parser.data["fchk_note"] = f"Orbitals from {name}."
    return True


def _add_sibling_output(parser, fchk_path, progress):
    """Opening a checkpoint also reads the job's output file when it is beside it."""
    log_path = _first_existing(fchk_path, (".log", ".out"))
    if not log_path:
        return
    try:
        text = read_gaussian_text(log_path)
    except OSError as exc:
        logging.warning("Reading the output file %s failed: %s", log_path, exc)
        return
    if "Entering Gaussian System" not in text[:4096] and (
        "Gaussian, Inc." not in text[:4096]
    ):
        return
    combined = GaussianParser()
    combined.load_from_memory(text, log_path, progress=progress)
    if attach_fchk(combined, fchk_path):
        parser.filename = combined.filename
        parser.raw_content = combined.raw_content
        parser.lines = combined.lines
        parser.data = combined.data
    else:
        parser.data["fchk_note"] = combined.data["fchk_note"]


def load_gaussian_parser(path, parent):
    """Read and parse *path* behind a progress dialog.

    Returns a loaded :class:`GaussianParser`, or ``None`` if the user cancelled.
    I/O errors propagate to the caller, which owns the error reporting.
    """
    progress = LoadProgress(parent, path)
    try:
        content = read_gaussian_text(path)
        if progress.cancelled:
            return None
        progress.read_done(path)

        parser = GaussianParser()
        try:
            if is_fchk_path(path) and looks_like_fchk(content):
                parser.load_fchk(content, path, progress=progress.on_parse_step)
                parser.data["fchk_path"] = path
                _add_sibling_output(parser, path, progress.on_parse_step)
            else:
                parser.load_from_memory(content, path, progress=progress.on_parse_step)
                fchk_path = sibling_fchk(path)
                if fchk_path:
                    attach_fchk(parser, fchk_path)
        except ParseCancelled:
            return None
        return parser
    finally:
        progress.close()
