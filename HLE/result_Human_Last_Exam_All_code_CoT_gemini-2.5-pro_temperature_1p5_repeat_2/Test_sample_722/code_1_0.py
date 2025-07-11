def analyze_sec_mals_data():
    """
    Analyzes the SEC-MALS experimental data step-by-step to determine the correct conclusion.
    """
    # Theoretical masses (in kDa)
    mass_A_theo = 25
    mass_B_theo = 150
    mass_C_theo = 60
    mass_D_theo = 100
    mass_kinase = 40

    print("--- Step-by-Step Analysis ---")

    # --- Experiment 1: Individual Proteins ---
    print("\n--- Analysis of Experiment 1: Individual Proteins ---")
    mass_A_exp1 = 50
    mass_B_exp1 = 300
    mass_C_exp1 = 60
    mass_D_exp1 = 100

    # Calculate oligomeric states
    print(f"Protein A: Theoretical mass is {mass_A_theo} kDa, experimental is {mass_A_exp1} kDa.")
    print(f"This means Protein A is a homodimer ({mass_A_theo} + {mass_A_theo} = {mass_A_exp1} kDa).")
    a_oligomer_mass = mass_A_exp1

    print(f"Protein B: Theoretical mass is {mass_B_theo} kDa, experimental is {mass_B_exp1} kDa.")
    print(f"This means Protein B is a homodimer ({mass_B_theo} + {mass_B_theo} = {mass_B_exp1} kDa).")
    b_oligomer_mass = mass_B_exp1

    print(f"Protein C: Theoretical mass is {mass_C_theo} kDa, experimental is {mass_C_exp1} kDa. It is a monomer.")
    c_oligomer_mass = mass_C_exp1
    
    print(f"Protein D: Theoretical mass is {mass_D_theo} kDa, experimental is {mass_D_exp1} kDa. It is a monomer.")
    d_oligomer_mass = mass_D_exp1

    # --- Experiment 2: Mixture of A, B, C, D ---
    print("\n--- Analysis of Experiment 2: All four proteins mixed ---")
    peak1_exp2 = 300
    peak2_exp2 = 210
    print(f"Detected peaks are at {peak1_exp2} kDa and {peak2_exp2} kDa.")
    print(f"The {peak1_exp2} kDa peak corresponds to the free Protein B homodimer ({b_oligomer_mass} kDa).")
    
    # Calculate mass of the other complex
    complex_ACD_mass = a_oligomer_mass + c_oligomer_mass + d_oligomer_mass
    print(f"The second peak at {peak2_exp2} kDa corresponds to a complex of Protein A dimer + Protein C + Protein D.")
    print(f"Calculation: {a_oligomer_mass} (A-dimer) + {c_oligomer_mass} (C) + {d_oligomer_mass} (D) = {complex_ACD_mass} kDa.")
    print("Conclusion: Non-phosphorylated Protein A (as a dimer) has a higher affinity for C+D than the Protein B dimer does.")

    # --- Experiment 3: Mixture with Kinase ---
    print("\n--- Analysis of Experiment 3: All four proteins mixed with Kinase ---")
    peak1_exp3 = 25
    peak2_exp3 = 40
    peak3_exp3 = 460
    print(f"Detected peaks are at {peak1_exp3} kDa, {peak2_exp3} kDa, and {peak3_exp3} kDa.")
    print(f"The {peak1_exp3} kDa peak is Protein A monomer ({mass_A_theo} kDa). Phosphorylation by the kinase broke the A-dimer.")
    print(f"The {peak2_exp3} kDa peak is the free kinase ({mass_kinase} kDa).")

    # Calculate mass of the large complex
    complex_BCD_mass = b_oligomer_mass + c_oligomer_mass + d_oligomer_mass
    print(f"The large peak at {peak3_exp3} kDa corresponds to a complex of Protein B dimer + Protein C + Protein D.")
    print(f"Calculation: {b_oligomer_mass} (B-dimer) + {c_oligomer_mass} (C) + {d_oligomer_mass} (D) = {complex_BCD_mass} kDa.")
    print("Conclusion: Phosphorylation of Protein A decreases its affinity for C+D, allowing the Protein B dimer to form a complex with C+D instead.")
    
    # --- Experiment 4: Dephosphorylation ---
    print("\n--- Analysis of Experiment 4: Sample from Exp 3 after dephosphorylation ---")
    peak1_exp4 = 50
    peak2_exp4 = 460
    print(f"Detected peaks are {peak1_exp4} kDa and {peak2_exp4} kDa.")
    print(f"The {peak1_exp4} kDa peak corresponds to the reformed Protein A dimer ({a_oligomer_mass} kDa).")
    print(f"The {peak2_exp4} kDa peak is the stable complex of B-dimer + C + D ({complex_BCD_mass} kDa) that formed in Experiment 3.")
    print("Conclusion: The complex of B+C+D, once formed, is very stable. Even when Protein A is dephosphorylated, it cannot displace the B-dimer from the complex. This supports the conclusion from Exp 2 that non-phosphorylated A has higher affinity in a competitive scenario.")

    # --- Final Evaluation ---
    print("\n--- Final Evaluation of Answer Choices ---")
    print("Based on the analysis:")
    print("1. From Exp 2, non-phosphorylated Protein A has a higher binding affinity for C+D than Protein B does.")
    print("2. From Exp 1, 2, 3 and 4, Protein B is consistently observed as a 300 kDa dimer or as part of a complex derived from the dimer.")
    print("3. Protein A can exist as a monomer (Exp 3, 25 kDa) or a dimer (Exp 1, 2, 4, 50 kDa). Proteins C and D are monomers.")
    print("\nThis matches statement G: 'Protein B never has a higher affinity for protein C or D than nonphosphorylated protein A, protein B always exists as a dimer in opposition to proteins A, C, and D.'")

analyze_sec_mals_data()