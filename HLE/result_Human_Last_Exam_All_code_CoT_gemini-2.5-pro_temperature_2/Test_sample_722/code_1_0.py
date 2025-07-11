def solve_protein_puzzle():
    """
    Analyzes the SEC-MALS experiments to determine protein interactions.
    """

    # --- Protein Masses ---
    # Theoretical (monomer) mass from sequence
    mass_A_mono = 25
    mass_B_mono = 150
    mass_C_mono = 60
    mass_D_mono = 100
    mass_kinase = 40

    # --- Step 1: Analyze Experiment 1 (Individual Proteins) ---
    print("--- Analysis of Experiment 1: Oligomeric States ---")
    mass_A_exp1 = 50
    mass_B_exp1 = 300
    mass_C_exp1 = 60
    mass_D_exp1 = 100
    
    state_A = "Dimer" if mass_A_exp1 / mass_A_mono == 2 else "Monomer"
    state_B = "Dimer" if mass_B_exp1 / mass_B_mono == 2 else "Monomer"
    state_C = "Monomer" # mass_C_exp1 == mass_C_mono
    state_D = "Monomer" # mass_D_exp1 == mass_D_mono

    print(f"Protein A (Theoretical: {mass_A_mono} kDa) is a {state_A} ({mass_A_exp1} kDa)")
    print(f"Protein B (Theoretical: {mass_B_mono} kDa) is a {state_B} ({mass_B_exp1} kDa)")
    print(f"Protein C (Theoretical: {mass_C_mono} kDa) is a {state_C} ({mass_C_exp1} kDa)")
    print(f"Protein D (Theoretical: {mass_D_mono} kDa) is a {state_D} ({mass_D_exp1} kDa)")
    print("\n")
    
    # --- Step 2: Analyze Experiment 2 (A, B, C, D Mixed) ---
    print("--- Analysis of Experiment 2: Unphosphorylated Mixture ---")
    peak1_exp2 = 300
    peak2_exp2 = 210
    
    complex_ACD = mass_A_exp1 + mass_C_exp1 + mass_D_exp1
    
    print(f"Detected peaks at {peak1_exp2} kDa and {peak2_exp2} kDa.")
    print(f"The {peak1_exp2} kDa peak corresponds to the free Protein B dimer.")
    print(f"The {peak2_exp2} kDa peak corresponds to the A-dimer, C, and D complex:")
    print(f"Equation: Mass(A-dimer) + Mass(C) + Mass(D) = {mass_A_exp1} + {mass_C_exp1} + {mass_D_exp1} = {complex_ACD} kDa.")
    print("Conclusion: Dephosphorylated Protein A has a higher affinity for C and/or D than Protein B does.\n")
    
    # --- Step 3: Analyze Experiment 3 (Mixed with Kinase) ---
    print("--- Analysis of Experiment 3: Mixture with Kinase (Phosphorylation) ---")
    peak1_exp3 = 25
    peak2_exp3 = 40
    peak3_exp3 = 460
    
    complex_BCD = mass_B_exp1 + mass_C_exp1 + mass_D_exp1
    
    print(f"Detected peaks at {peak1_exp3} kDa, {peak2_exp3} kDa, and {peak3_exp3} kDa.")
    print(f"The {peak1_exp3} kDa peak is Protein A as a monomer.")
    print(f"The {peak2_exp3} kDa peak is the free kinase.")
    print(f"The {peak3_exp3} kDa peak corresponds to the B-dimer, C, and D complex:")
    print(f"Equation: Mass(B-dimer) + Mass(C) + Mass(D) = {mass_B_exp1} + {mass_C_exp1} + {mass_D_exp1} = {complex_BCD} kDa.")
    print("Conclusion: Phosphorylation (of Protein A) causes it to dissociate from C/D and also break into monomers. This allows Protein B to bind C and D.\n")
    
    # --- Step 4: Final Evaluation ---
    print("--- Final Evaluation and Answer Selection ---")
    print("Based on the analysis:")
    print("1. Exp 2 shows that nonphosphorylated Protein A has a higher affinity for C/D than Protein B does.")
    print("2. All experiments show that Protein B is consistently a 300 kDa dimer.")
    print("3. Protein B's state (always dimer) is different from Protein A (dimer or monomer), Protein C (monomer), and Protein D (monomer).")
    print("This directly supports all clauses of statement G.")
    print("\nMost Correct Answer: G. Protein B never has a higher affinity for protein C or D than nonphosphorylated protein A, protein B always exists as a dimer in opposition to proteins A, C, and D.")

solve_protein_puzzle()
<<<G>>>