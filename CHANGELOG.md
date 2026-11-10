# Changelog

## 0.2.0 — 2026-11-11

- Read every orbital energy of a job that prints no MO table. Only a few were read, and with wrong values.
- Read thermochemistry from `#T` logs, and the total energy of a double hybrid such as B2PLYP.
- Do not report the `<S**2>` of a stability analysis as spin contamination of a singlet.
- Point the forces of a molecule above 50 atoms the right way; Gaussian prints no input orientation there.
- Remove the **Load .fchk** button: the one **Open** button opens logs and formatted checkpoints, and a checkpoint beside a log is picked up by itself.
- The file box uses the same blue as the ORCA Result Analyzer.
- The description no longer says under development.

## 0.1.3 — 2026-11-10

- On MoleditPy 5.2.0 and later, pass another program's file on to the next plugin by returning `False` instead of raising, so the host no longer logs a traceback for it. Older hosts keep the exception.

## 0.1.2 — 2026-11-08

- Show every optimization cycle by default when opening the optimization tab.
- Keep the scan-point view as the default for relaxed scans with multiple points.
- Add regressions using real Gaussian logs, formatted checkpoints and real Qt windows.
