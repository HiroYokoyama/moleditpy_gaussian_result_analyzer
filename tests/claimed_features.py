"""The features the Gaussian analyzer claims, with the real output that shows each."""

# Every claimed window: (fixture, main-window slot, window class).
CLAIMS = [
    ("water_6d_mo", "show_scf_trace", "SCFTraceDialog"),
    ("water_6d_mo", "show_mo_analyzer", "MODialog"),
    ("h2co_5d7f_mo", "show_mo_analyzer", "MODialog"),
    ("o2_triplet_uhf_opt", "show_mo_analyzer", "MODialog"),
    ("o2_triplet_uhf_opt", "show_trajectory", "TrajectoryResultDialog"),
    ("ethane_scan", "show_trajectory", "TrajectoryResultDialog"),
    ("benzene_gradient", "show_forces", "ForceViewerDialog"),
    ("water_chelpg", "show_charges", "ChargeDialog"),
    ("water_mk", "show_charges", "ChargeDialog"),
    ("water_6d_mo", "show_dipole", "DipoleDialog"),
    ("water_raman", "show_freq", "FrequencyDialog"),
    ("benzene_opt_freq", "show_thermal", "ThermalTableDialog"),
    ("h2co_tddft", "show_tddft", "TDDFTDialog"),
    ("methanol_nmr", "show_nmr", "NMRDialog"),
    ("water_6d_mo", "show_bond_analysis", "BondAnalysisDialog"),
    ("water_mp2", "show_energy_components", "EnergyComponentsDialog"),
    ("water_ccsdt", "show_energy_components", "EnergyComponentsDialog"),
    ("water_d3bj", "show_properties", "PropertiesDialog"),
    ("water_6d_mo", "show_properties", "PropertiesDialog"),
]

# (fixture, slot, controls that must be there, controls that must not).
CONTROLS = [
    ("ethane_scan", "show_trajectory", {"Scan Points", "All Trajectory", "Coordinate"}, {"Load MEP TRJ"}),
    ("h2co_tddft", "show_tddft", {"Absorption", "CD"}, set()),
    ("methanol_nmr", "show_nmr", {"Simulate Coupling"}, set()),
    ("water_raman", "show_freq", {"Open IR/Raman Spectrum..."}, {"Show Dipole Change"}),
]

TABS = [
    ("water_6d_mo", "show_bond_analysis", "NBO Orbitals", "Mayer"),
    ("water_raman", "show_freq", "Mode Vectors", "Dipole"),
]

# (fixture, button attribute): a job that has the result enables the button.
BUTTONS = [
    ("water_6d_mo", "btn_scf"),
    ("water_6d_mo", "btn_mo"),
    ("o2_triplet_uhf_opt", "btn_traj"),
    ("ethane_scan", "btn_traj"),
    ("benzene_gradient", "btn_forces"),
    ("water_chelpg", "btn_charge"),
    ("water_6d_mo", "btn_dipole"),
    ("water_raman", "btn_freq"),
    ("benzene_opt_freq", "btn_therm"),
    ("h2co_tddft", "btn_tddft"),
    ("methanol_nmr", "btn_nmr"),
]

# (fixture, button attribute): a job without the result leaves it disabled.
DISABLED = [
    ("water_6d_mo", "btn_freq"),
    ("water_6d_mo", "btn_tddft"),
    ("water_6d_mo", "btn_nmr"),
    ("water_6d_mo", "btn_traj"),
    ("water_6d_mo", "btn_therm"),
]
