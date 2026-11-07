"""Strip identifying and non-scientific text from a Gaussian log before committing it.

A log records where and by whom it was run: the account name and machine in the
archive block, and the full job path in the first lines. Nothing in the parser
reads these, so a fixture carries placeholders instead. The licence banner and
all timings are dropped too, and tests/test_fixture_privacy.py holds every
committed fixture to this.

Usage: python tools/sanitize_gaussian_log.py INPUT OUTPUT
"""

import re
import sys

_WINDOWS_PATH = re.compile(r"[A-Za-z]:\\[^\s\"]*")
_UNIX_HOME = re.compile(r"/(?:home|Users|scratch|work|tmp)/[^\s\"]*")
_SCRATCH_NAME = re.compile(r"Gau-\d+")
_ARCHIVE_START = re.compile(r"^ 1[\\|]1[\\|]")
_DROPPED_LINE = re.compile(
    r"^\s*(?:Job cpu time|Elapsed time|File lengths|Leave Link|Entering Link|"
    r"Input=|Output=|Initial command|\(Enter |Default route)"
)
_ARCHIVE_WIDTH = 70


def _clean_archive(block):
    """Replace the machine and user fields of one archive entry, keeping the rest."""
    text = "".join(line[1:] for line in block)
    separator = "\\" if text[3:4] == "\\" else "|"
    fields = text.split(separator)
    if len(fields) > 8:
        fields[2] = "HOST"
        fields[7] = "USER"
        fields[8] = "DATE"
    text = separator.join(fields)
    return [
        " " + text[i : i + _ARCHIVE_WIDTH] for i in range(0, len(text), _ARCHIVE_WIDTH)
    ]


def sanitize(text):
    """Return *text* with identifying paths, names, banner and timings removed."""
    lines = text.splitlines()
    out = []
    i = 0
    # Everything before the first "****" rule of the route echo is banner.
    while i < len(lines) and not lines[i].startswith(" Gaussian "):
        i += 1
    if i < len(lines):
        out.append(lines[i])
        i += 1
    while i < len(lines):
        line = lines[i]
        if _ARCHIVE_START.match(line):
            block = []
            while i < len(lines) and lines[i].strip():
                block.append(lines[i])
                if lines[i].rstrip().endswith("@"):
                    i += 1
                    break
                i += 1
            out.extend(_clean_archive(block))
            # What follows the archive is the closing quotation, up to the times.
            while i < len(lines) and not lines[i].lstrip().startswith(
                ("Job cpu time", "Normal termination", "Error termination")
            ):
                i += 1
            continue
        i += 1
        if _DROPPED_LINE.match(line):
            continue
        line = _WINDOWS_PATH.sub("PATH", line)
        line = _UNIX_HOME.sub("PATH", line)
        line = _SCRATCH_NAME.sub("Gau-0", line)
        out.append(line)
    return "\n".join(out) + "\n"


def main(argv):
    """Sanitize argv[1] into argv[2]."""
    with open(argv[1], encoding="utf-8", errors="replace") as handle:
        text = handle.read()
    with open(argv[2], "w", encoding="utf-8", newline="\n") as handle:
        handle.write(sanitize(text))


if __name__ == "__main__":
    main(sys.argv)
