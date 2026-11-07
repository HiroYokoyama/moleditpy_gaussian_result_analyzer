"""NMR dialog mixin: exporting the shielding/coupling tables and spectrum plot."""

import os
import csv
from PySide6.QtWidgets import QApplication, QMessageBox, QFileDialog
from .utils import get_default_export_path, notify


class _NMRExportMixin:
    def export_spectrum(self):
        """Export spectrum to file"""
        current_nucleus = self.current_nucleus
        default_path = get_default_export_path(
            self.file_path, suffix=f"_nmr_{current_nucleus}_spectrum", extension=".png"
        )

        filename, _ = QFileDialog.getSaveFileName(
            self,
            "Export Spectrum",
            default_path,
            "PNG Image (*.png);;PDF (*.pdf);;SVG (*.svg)",
        )

        if filename:
            try:
                self.figure.savefig(filename, dpi=300, bbox_inches="tight")
                notify(
                    self, f"Spectrum exported to: {os.path.basename(filename)}", 5000
                )
            # ValueError: savefig rejects an unsupported extension typed in
            # the file dialog, which previously escaped the slot.
            except (OSError, ValueError) as e:
                QMessageBox.critical(self, "Error", f"Export failed:\n{e}")

    def export_spectrum_csv(self):
        """Export spectrum data to CSV (Stick or Real)"""
        if not self.displayed_data:
            QMessageBox.warning(self, "No Data", "No spectrum data to export.")
            return

        default_path = get_default_export_path(
            self.file_path,
            suffix=f"_nmr_{self.current_nucleus}_sticks",
            extension=".csv",
        )
        filename, _ = QFileDialog.getSaveFileName(
            self, "Export Spectrum CSV", default_path, "CSV Files (*.csv)"
        )
        if not filename:
            return

        try:
            with open(filename, "w", encoding="utf-8") as f:
                # Check mode
                is_real = (
                    getattr(self, "chk_real_spectrum", None) is not None
                    and self.chk_real_spectrum.isChecked()
                )

                if is_real:
                    # Export XY data from matplotlib line
                    ax = self.figure.axes[0] if self.figure.axes else None
                    if ax and ax.lines:
                        # The last line drawn is usually the spectrum curve or baseline?
                        # In plot_real_spectrum, we draw: ax.plot(x_hz_grid / spectrometer_freq, y_total, 'b-', linewidth=1.2)
                        # We should make sure we grab the right line.
                        # It is the only 'b-' line usually.

                        # Better approach: Re-calculate the data?
                        # Extracting from plot is hacky but consistent with what is seen.
                        # Let's try to extract x,y data from the first non-stem line.

                        line = None
                        for ln in ax.lines:
                            # Stems are Line2D but usually handled differently.
                            # nmrsim plot is a standard plot.
                            if ln.get_marker() == "None" and ln.get_linestyle() == "-":
                                line = ln
                                # If we have multiple, the 'blue' one is the spectrum.
                                if ln.get_color() == "b":
                                    break

                        if line:
                            x_data = line.get_xdata()
                            y_data = line.get_ydata()

                            f.write("Chemical Shift (ppm),Intensity\n")
                            for x, y in zip(x_data, y_data):
                                f.write(f"{x:.6f},{y:.6f}\n")
                        else:
                            f.write("Error: Could not extract curve data.\n")
                    else:
                        f.write("Error: No plot found.\n")
                else:
                    # Stick Spectrum Export
                    # Export Peak List: Shift, Intensity
                    f.write("Chemical Shift (ppm),Intensity,AtomIndices\n")

                    # We can use peaks_metadata if available, or reconstruct
                    if (
                        getattr(self, "peaks_metadata", None) is not None
                        and self.peaks_metadata
                    ):
                        for (
                            shift,
                            intensity,
                            _is_merged,
                            atom_indices,
                        ) in self.peaks_metadata:
                            indices_str = ";".join(str(i) for i in atom_indices)
                            f.write(f"{shift:.6f},{intensity:.4f},{indices_str}\n")
                    else:
                        # Fallback
                        f.write("# No peak data available.\n")

            notify(
                self, f"Spectrum data exported to: {os.path.basename(filename)}", 5000
            )

        # Qt slot: a slot must never crash the app (CONTRIBUTING.md 4B)
        except Exception as e:  # pylint: disable=broad-exception-caught
            QMessageBox.critical(self, "Error", f"Export failed:\n{e}")

    def export_table_csv(self):
        """Export table data to CSV"""
        current_nucleus = self.current_nucleus
        default_path = get_default_export_path(
            self.file_path, suffix=f"_nmr_{current_nucleus}_table", extension=".csv"
        )

        filename, _ = QFileDialog.getSaveFileName(
            self, "Export Table CSV", default_path, "CSV Files (*.csv)"
        )
        if not filename:
            return

        try:
            with open(filename, "w", newline="", encoding="utf-8") as f:
                # Headers
                headers = []
                for c in range(self.table.columnCount()):
                    headers.append(self.table.horizontalHeaderItem(c).text())
                writer = csv.writer(f)
                writer.writerow(headers)

                # Rows
                for r in range(self.table.rowCount()):
                    cols = []
                    for c in range(self.table.columnCount()):
                        it = self.table.item(r, c)
                        text = it.text() if it else ""
                        cols.append(text)
                    writer.writerow(cols)

            notify(self, f"Table data exported to: {os.path.basename(filename)}", 5000)
        # Qt slot: a slot must never crash the app (CONTRIBUTING.md 4B)
        except Exception as e:  # pylint: disable=broad-exception-caught
            QMessageBox.critical(self, "Error", f"Export failed:\n{e}")

    def get_j_coupling_string(self, atom_indices):
        """Format J-couplings for a list of atoms (merged or single)"""
        if not self.couplings:
            return ""

        # If merged group, we might have too many couplings.
        # Strategy:
        # 1. Collect all couplings involving ANY atom in atom_indices.
        # 2. Filter out couplings WITHIN the group (intra-group couplings often effectively 0 or infinite depending on equivalence,
        #    but usually we care about couplings to OUTSIDE).
        # 3. Format as "AtomSymIdx=J"
        # 4. If multiple atoms in group have SAME coupling to SAME partner (magnetic equivalence), show once.

        relevant_couplings = []  # (PartnerIdx, J)
        group_set = set(atom_indices)

        # Optimize: Pre-filter couplings? dataset is small usually.
        for c in self.couplings:
            idx1 = c["atom_idx1"]
            idx2 = c["atom_idx2"]
            j_val = abs(c["coupling"])

            if j_val < 0.1:
                continue  # Ignore tiny couplings

            partner = None
            if idx1 in group_set and idx2 not in group_set:
                partner = idx2
            elif idx2 in group_set and idx1 not in group_set:
                partner = idx1

            if partner is not None:
                relevant_couplings.append((partner, j_val))

        if not relevant_couplings:
            return ""

        # If merged group, we might see duplicates: (PartnerA, 7.5) from Atom1, (PartnerA, 7.5) from Atom2...
        # We should average them or check consistency.
        # Group by partner
        partner_map = {}
        for p, j in relevant_couplings:
            if p not in partner_map:
                partner_map[p] = []
            partner_map[p].append(j)

        # Format strings
        parts = []

        # Sort partners by Atom Index for consistency
        sorted_partners = sorted(partner_map.keys())

        for p in sorted_partners:
            # Find symbol for partner
            p_item = next((d for d in self.data if d.get("atom_idx", None) == p), None)
            p_sym = p_item.get("atom_sym", "") if p_item else ""
            p_label = f"{p_sym}{p}"

            j_list = partner_map[p]
            avg_j = sum(j_list) / len(j_list)

            # If standard deviation is high, maybe indicate?
            # For now, just show average.
            parts.append(f"{p_label}={avg_j:.1f}")

        return ", ".join(parts)

    def copy_table(self):
        """Copy table data to clipboard"""
        text = "Idx\tNucleus\tShielding\tShift\tJ-coupling\n"
        for r in range(self.table.rowCount()):
            cols = []
            for c in range(self.table.columnCount()):
                it = self.table.item(r, c)
                cols.append(it.text() if it else "")
            text += "\t".join(cols) + "\n"
        QApplication.clipboard().setText(text)
        notify(self, "Table data copied to clipboard!", 5000)
