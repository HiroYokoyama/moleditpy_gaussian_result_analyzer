# Changelog

## Unreleased

- On MoleditPy 5.2.0 and later, pass another program's file on to the next plugin by returning `False` instead of raising, so the host no longer logs a traceback for it. Older hosts keep the exception.

## 0.1.2 — 2026-11-08

- Show every optimization cycle by default when opening the optimization tab.
- Keep the scan-point view as the default for relaxed scans with multiple points.
- Add regressions using real Gaussian logs, formatted checkpoints and real Qt windows.
