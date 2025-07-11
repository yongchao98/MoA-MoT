def simulate_flow_cytometry_competition_assay():
    """
    This script simulates and explains the data from the MUC1 experiment
    to highlight the importance of the anti-FLAG antibody control.
    """

    print("--- Simulating Flow Cytometry Data ---")
    print("The experiment measures Mean Fluorescence Intensity (MFI) for two antibodies under two conditions.\n")

    # --- Simulated MFI Values ---

    # Condition 1: PBS (Control)
    # The anti-FLAG antibody confirms MUC1 is present on the cell surface.
    mfi_pbs_flag = 1200
    # The anti-glycoMUC1 antibody shows strong binding to its target.
    mfi_pbs_glycomucin = 950

    # Condition 2: 500 mM GalNAc (Inhibitor)
    # The anti-FLAG MFI should remain high, confirming MUC1 is still on the surface.
    mfi_galnac_flag = 1180
    # The anti-glycoMUC1 MFI should be very low, showing its binding was inhibited.
    mfi_galnac_glycomucin = 85

    # --- Analysis & Interpretation ---

    print(f"--- Results for Anti-FLAG Antibody (Control for MUC1 expression) ---")
    print(f"MFI in PBS: {mfi_pbs_flag}")
    print(f"MFI in 500 mM GalNAc: {mfi_galnac_flag}")
    expression_ratio = mfi_galnac_flag / mfi_pbs_flag
    print(f"Analysis: The amount of FLAG-MUC1 on the surface in the GalNAc condition is {expression_ratio:.2%} of the PBS control ({mfi_galnac_flag} / {mfi_pbs_flag}).")
    print("Conclusion: The high concentration of GalNAc did not significantly alter the surface expression of the MUC1 protein.\n")


    print(f"--- Results for Anti-glycoMUC1 Antibody (Antibody of interest) ---")
    print(f"MFI in PBS: {mfi_pbs_glycomucin}")
    print(f"MFI in 500 mM GalNAc: {mfi_galnac_glycomucin}")
    inhibition_percentage = (1 - (mfi_galnac_glycomucin / mfi_pbs_glycomucin)) * 100
    print(f"Analysis: The antibody binding was inhibited by {inhibition_percentage:.2f}%. Calculation: (1 - ({mfi_galnac_glycomucin} / {mfi_pbs_glycomucin})) * 100")
    print("Conclusion: Free GalNAc effectively competed with the MUC1 glycan for antibody binding.\n")

    print("--- Final Justification ---")
    print("Because the anti-FLAG control confirmed that MUC1 surface levels were stable, we can be confident that the dramatic drop in the anti-glycoMUC1 signal is due to specific inhibition, not a loss of the target protein.")
    print("This demonstrates why it is essential to include the anti-FLAG antibody with the primary antibodies: to verify that GalNAc has not altered the surface expression of MUC1.")


if __name__ == "__main__":
    simulate_flow_cytometry_competition_assay()