"""Every committed fixture must be free of the machine, the account and the paths."""

import os
import re
import unittest

SAMPLES = os.path.join(os.path.dirname(__file__), "sample_outputs")

_FORBIDDEN = {
    "Windows path": re.compile(r"[A-Za-z]:\\[A-Za-z0-9_.]"),
    "Unix home": re.compile(r"/(?:home|Users|scratch)/"),
    "scratch file": re.compile(r"Gau-[1-9]"),
    "job path echo": re.compile(r"^\s*(?:Input=|Output=|Initial command)", re.M),
    "timing": re.compile(r"Job cpu time|Elapsed time"),
    "UUID": re.compile(r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}"),
    "private address": re.compile(r"\b(?:192\.168|10\.\d+\.\d+)\.\d+\b"),
}
_ARCHIVE = re.compile(r"^ 1[\\|]1[\\|]")


def _samples():
    return sorted(f for f in os.listdir(SAMPLES) if f.endswith((".log", ".out")))


_SAFE_MOLECULES = {"benzene", "ethane", "h2co", "heptadecane", "methyl", "o2", "water"}


def _checkpoints():
    return sorted(f for f in os.listdir(SAMPLES) if f.endswith(".fchk"))


def _read(name):
    with open(os.path.join(SAMPLES, name), encoding="utf-8") as handle:
        return handle.read()


class TestFixturePrivacy(unittest.TestCase):
    def test_there_are_fixtures_to_check(self):
        self.assertGreater(len(_samples()), 10)

    def test_no_identifying_text(self):
        for name in _samples():
            text = _read(name)
            for label, pattern in _FORBIDDEN.items():
                with self.subTest(fixture=name, kind=label):
                    self.assertIsNone(pattern.search(text))

    def test_checkpoints_hold_no_identifying_text(self):
        self.assertGreater(len(_checkpoints()), 10)
        for name in _checkpoints():
            text = _read(name)
            self.assertIn(name.split("_")[0], _SAFE_MOLECULES)
            for label, pattern in _FORBIDDEN.items():
                with self.subTest(fixture=name, kind=label):
                    self.assertIsNone(pattern.search(text))

    def test_archive_entries_carry_placeholders(self):
        for name in _samples():
            lines = _read(name).splitlines()
            for i, line in enumerate(lines):
                if not _ARCHIVE.match(line):
                    continue
                with self.subTest(fixture=name):
                    head = "".join(row[1:] for row in lines[i : i + 3])
                    fields = re.split(r"[\\|]", head)
                    self.assertEqual(fields[2], "HOST")
                    self.assertEqual(fields[7], "USER")
                    self.assertEqual(fields[8], "DATE")


if __name__ == "__main__":
    unittest.main()
