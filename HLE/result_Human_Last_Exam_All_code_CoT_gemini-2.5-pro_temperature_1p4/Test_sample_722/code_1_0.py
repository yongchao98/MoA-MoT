def analyze_protein_interactions():
    """
    Analyzes SEC-MALS data to determine protein interactions and evaluates conclusions.
    """

    # --- Masses ---
    mass_A_theo = 25
    mass_B_theo = 150
    mass_C_theo = 60
    mass_D_theo = 100
    mass_kinase = 40

    print("--- Step 1: Analysis of Experiment 1 (Individual Proteins) ---")
    mass_A_exp1 = 50
    mass_B_exp1 = 300
    mass_C_exp1 = 60
    mass_D_exp1 = 100
    
    oligomer_A = mass_A_exp1 / mass_A_theo
    oligomer_B = mass_B_exp1 / mass_B_theo
    
    print(f"Protein A: Experimental mass {mass_A_exp1} kDa / Theoretical mass {mass_A_theo} kDa = {oligomer_A}. Conclusion: Protein A is a dimer (A2).")
    print(f"Protein B: Experimental mass {mass_B_exp1} kDa / Theoretical mass {mass_B_theo} kDa = {oligomer_B}. Conclusion: Protein B is a dimer (B2).")
    print(f"Protein C: Experimental mass {mass_C_exp1} kDa. Conclusion: Protein C is a monomer.")
    print(f"Protein D: Experimental mass {mass_D_exp1} kDa. Conclusion: Protein D is a monomer.\n")

    print("--- Step 2: Analysis of Experiment 2 (A+B+C+D Mixed) ---")
    peak1_exp2 = 300
    peak2_exp2 = 210
    
    print(f"Observed peak at {peak1_exp2} kDa corresponds to the Protein B dimer (B2), which is unbound.")
    
    complex_ACD_mass = mass_A_exp1 + mass_C_exp1 + mass_D_exp1
    print(f"Let's test if the second peak ({peak2_exp2} kDa) is a complex of A2, C, and D:")
    print(f"Mass(A2) + Mass(C) + Mass(D) = {mass_A_exp1} + {mass_C_exp1} + {mass_D_exp1} = {complex_ACD_mass} kDa.")
    print("This matches perfectly. Conclusion: A2 binds to C and D, outcompeting B2. So, Affinity(A2) > Affinity(B2) for C and D.\n")

    print("--- Step 3: Analysis of Experiment 3 (A+B+C+D + Kinase) ---")
    peak1_exp3 = 25
    peak2_exp3 = 40
    peak3_exp3 = 460
    
    print(f"Observed peak at {peak1_exp3} kDa matches the theoretical mass of Protein A monomer. This implies phosphorylation by the kinase breaks the A2 dimer.")
    print(f"Observed peak at {peak2_exp3} kDa is the free kinase.")
    
    complex_BCD_mass = mass_B_exp1 + mass_C_exp1 + mass_D_exp1
    print(f"Let's test if the third peak ({peak3_exp3} kDa) is a complex of B2, C, and D:")
    print(f"Mass(B2) + Mass(C) + Mass(D) = {mass_B_exp1} + {mass_C_exp1} + {mass_D_exp1} = {complex_BCD_mass} kDa.")
    print("This matches perfectly. Conclusion: Phosphorylation of A decreases its affinity for C and D, allowing B2 to bind them.\n")

    print("--- Step 4: Analysis of Experiment 4 (Dephosphorylation) ---")
    peak1_exp4 = 50
    peak2_exp4 = 460
    print(f"Observed peak at {peak1_exp4} kDa shows that dephosphorylated Protein A has re-formed its dimer (A2).")
    print(f"Observed peak at {peak2_exp4} kDa shows the B2CD complex is stable and is not broken up by A2.\n")
    
    print("--- Step 5: Final Conclusion ---")
    print("Based on the analysis:")
    print("1. From Exp 2, nonphosphorylated Protein A has a higher affinity for C and D than Protein B.")
    print("2. From all experiments, Protein B is always a dimer, while Protein A can be a monomer or dimer, and C and D are monomers.")
    print("These two points directly support answer choice G.\n")
    
# Run the analysis
analyze_protein_interactions()

# Final Answer
print("<<<G>>>")