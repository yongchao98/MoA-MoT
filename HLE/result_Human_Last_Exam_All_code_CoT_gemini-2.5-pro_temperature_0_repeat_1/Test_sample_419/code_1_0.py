def simulate_binding_experiment():
    """
    This script simulates the results of a cell surface binding experiment
    to illustrate the importance of the anti-FLAG antibody control.
    """
    # The concentration of the inhibitor used in the experiment.
    inhibitor_concentration_mM = 500

    print("### Analyzing the Role of the Anti-FLAG Control ###\n")
    print(f"The experiment tests if {inhibitor_concentration_mM} mM GalNAc inhibits antibody binding through competition.\n")

    # --- SCENARIO 1: Ideal Result (Validating Competition) ---
    print("--- SCENARIO 1: The control confirms the hypothesis ---")
    # Hypothetical fluorescence signals (Arbitrary Fluorescence Units, AFU)
    pbs_condition = {'anti_glycan_MUC1_signal': 950, 'anti_FLAG_signal': 1100}
    galnac_condition_good = {'anti_glycan_MUC1_signal': 120, 'anti_FLAG_signal': 1080}

    print(f"Condition: PBS Control")
    print(f"  - Anti-Glycan-MUC1 Signal: {pbs_condition['anti_glycan_MUC1_signal']} AFU")
    print(f"  - Anti-FLAG Signal: {pbs_condition['anti_FLAG_signal']} AFU\n")

    print(f"Condition: {inhibitor_concentration_mM} mM GalNAc")
    print(f"  - Anti-Glycan-MUC1 Signal: {galnac_condition_good['anti_glycan_MUC1_signal']} AFU (Signal decreased)")
    print(f"  - Anti-FLAG Signal: {galnac_condition_good['anti_FLAG_signal']} AFU (Signal is stable)\n")

    print("Conclusion for Scenario 1:")
    print("The Anti-FLAG signal is stable, proving that the total amount of MUC1 on the cell surface has NOT changed.")
    print("Therefore, the decrease in the Anti-Glycan-MUC1 signal is correctly attributed to binding competition.\n")
    print("="*60 + "\n")

    # --- SCENARIO 2: Confounded Result (Invalidating the Conclusion) ---
    print("--- SCENARIO 2: The control reveals a confounding variable ---")
    # Here, the high GalNAc concentration is toxic and reduces MUC1 surface expression
    galnac_condition_bad = {'anti_glycan_MUC1_signal': 110, 'anti_FLAG_signal': 150}

    print(f"Condition: PBS Control")
    print(f"  - Anti-Glycan-MUC1 Signal: {pbs_condition['anti_glycan_MUC1_signal']} AFU")
    print(f"  - Anti-FLAG Signal: {pbs_condition['anti_FLAG_signal']} AFU\n")

    print(f"Condition: {inhibitor_concentration_mM} mM GalNAc")
    print(f"  - Anti-Glycan-MUC1 Signal: {galnac_condition_bad['anti_glycan_MUC1_signal']} AFU (Signal decreased)")
    print(f"  - Anti-FLAG Signal: {galnac_condition_bad['anti_FLAG_signal']} AFU (Signal ALSO decreased!)\n")

    print("Conclusion for Scenario 2:")
    print("The Anti-FLAG signal has also dropped significantly. This reveals that the {inhibitor_concentration_mM} mM GalNAc treatment caused a loss of MUC1 protein from the cell surface.")
    print("The experiment is invalid because we cannot distinguish between competition and loss of target.\n")

    print("This simulation demonstrates why the anti-FLAG antibody must be added with the primary antibodies: to verify that GalNAc has not altered the surface expression of MUC1.")

# Run the simulation
simulate_binding_experiment()