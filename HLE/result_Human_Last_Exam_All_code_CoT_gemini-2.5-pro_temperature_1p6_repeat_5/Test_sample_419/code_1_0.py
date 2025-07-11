def analyze_experiment(anti_muc1_glycan_signal, anti_flag_signal, condition):
    """
    Analyzes the signal from a simulated flow cytometry or ELISA experiment.

    Args:
        anti_muc1_glycan_signal (float): The signal from the antibody that binds the TN-antigen on MUC1.
        anti_flag_signal (float): The signal from the anti-flag antibody, which measures MUC1 protein level.
        condition (str): The experimental condition, e.g., 'PBS Control' or 'GalNAc Treatment'.
    """
    print(f"--- Analyzing Condition: {condition} ---")
    print(f"Signal from Anti-MUC1-Glycan Antibody: {anti_muc1_glycan_signal:.1f}")
    print(f"Signal from Anti-FLAG Antibody (Surface MUC1 Control): {anti_flag_signal:.1f}\n")
    return anti_muc1_glycan_signal, anti_flag_signal

def main():
    """
    Simulates and interprets an antibody inhibition experiment.
    """
    print("Simulating an experiment to test antibody specificity to glycosylated MUC1.\n")

    # Experimental parameters from the problem description
    galnac_concentration = 500  # in mM

    # --- SIMULATED DATA ---
    # Case 1: Ideal result where GalNAc does not affect MUC1 expression.
    print("### SCENARIO 1: The experiment is controlled correctly. ###")
    # In the PBS control condition, both antibodies give a strong signal.
    pbs_muc1_signal, pbs_flag_signal = analyze_experiment(100.0, 100.0, "PBS Control")

    # In the GalNAc condition, the anti-MUC1-glycan signal drops, but the anti-flag signal is stable.
    galnac_treatment_str = f"{galnac_concentration} mM GalNAc Treatment"
    galnac_muc1_signal, galnac_flag_signal = analyze_experiment(15.0, 98.0, galnac_treatment_str)

    # Interpretation Logic
    # We check if the MUC1 protein level (measured by anti-flag) is stable.
    # Allow for a small (e.g., 10%) experimental variance.
    if galnac_flag_signal > (pbs_flag_signal * 0.9):
        print("INTERPRETATION:")
        print(f"The Anti-FLAG signal is stable, confirming that the total amount of MUC1 on the cell surface was NOT significantly altered by the {galnac_concentration} mM GalNAc treatment.")
        print("The decrease in the Anti-MUC1-Glycan signal can be confidently attributed to competitive inhibition by free GalNAc.")
        print("CONCLUSION: The experiment successfully shows the antibody is specific to the TN antigen.\n")
    else:
        print("INTERPRETATION:")
        print(f"The Anti-FLAG signal has decreased significantly, indicating that the {galnac_concentration} mM GalNAc treatment reduced MUC1 surface expression.")
        print("CONCLUSION: The experiment is inconclusive. The drop in Anti-MUC1-Glycan signal could be due to inhibition OR simply less protein on the surface.\n")


    # ---
    # Case 2: Hypothetical result where GalNAc is cytotoxic and reduces MUC1 expression.
    print("\n### SCENARIO 2: A hypothetical, poorly controlled result. ###")
    # PBS control is the same.
    pbs_muc1_signal_2, pbs_flag_signal_2 = analyze_experiment(100.0, 100.0, "PBS Control")
    # In this hypothetical case, the GalNAc treatment reduces BOTH signals.
    galnac_muc1_signal_2, galnac_flag_signal_2 = analyze_experiment(15.0, 20.0, galnac_treatment_str)
    
    if galnac_flag_signal_2 > (pbs_flag_signal_2 * 0.9):
         # This part will not be reached in scenario 2
        pass
    else:
        print("INTERPRETATION:")
        print(f"The Anti-FLAG signal has decreased significantly, indicating that the {galnac_concentration} mM GalNAc treatment likely reduced MUC1 surface expression.")
        print("CONCLUSION: The experiment is inconclusive. Without the anti-flag control, one might have wrongly concluded competitive inhibition was the only factor.")

if __name__ == "__main__":
    main()
