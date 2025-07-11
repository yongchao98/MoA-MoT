def solve_protein_puzzle():
    """
    Analyzes the SEC-MALS experiment data to determine the correct conclusion.
    The analysis is printed step-by-step.
    """
    
    # Protein theoretical masses
    mass_A_mono = 25
    mass_B_mono = 150
    mass_C_mono = 60
    mass_D_mono = 100
    mass_kinase = 40

    print("--- Step 1: Analysis of Experiment 1 (Individual Proteins) ---")
    print("This step determines the native oligomeric state of each protein.")
    
    # Protein A
    obs_mass_A = 50
    print(f"Protein A: Theoretical mass is {mass_A_mono} kDa, observed mass is {obs_mass_A} kDa.")
    print(f"Conclusion: Protein A exists as a homodimer (2 molecules of A).")
    print(f"Equation: {mass_A_mono} * 2 = {mass_A_mono * 2}")
    mass_A_native = obs_mass_A
    
    # Protein B
    obs_mass_B = 300
    print(f"\nProtein B: Theoretical mass is {mass_B_mono} kDa, observed mass is {obs_mass_B} kDa.")
    print(f"Conclusion: Protein B exists as a homodimer (2 molecules of B).")
    print(f"Equation: {mass_B_mono} * 2 = {mass_B_mono * 2}")
    mass_B_native = obs_mass_B

    # Protein C and D
    print(f"\nProtein C (60 kDa) and D (100 kDa) have observed masses equal to their theoretical masses, so they are monomers.")

    print("\n--- Step 2: Analysis of Experiment 2 (Proteins A, B, C, D Mixed) ---")
    print("This step determines binding preferences in the absence of modification.")
    
    peak1_exp2 = 300
    peak2_exp2 = 210
    
    print(f"Observed peaks are at {peak1_exp2} kDa and {peak2_exp2} kDa.")
    print(f"The {peak1_exp2} kDa peak matches the mass of the Protein B dimer ({mass_B_native} kDa). This means Protein B is not part of the other complex.")
    
    complex_ACD_mass = mass_A_native + mass_C_mono + mass_D_mono
    print(f"The {peak2_exp2} kDa peak corresponds to a complex of the remaining proteins: Protein A dimer + Protein C + Protein D.")
    print(f"Equation: {mass_A_native} + {mass_C_mono} + {mass_D_mono} = {complex_ACD_mass}")
    print("Conclusion: When competing, the Protein A dimer binds to Proteins C and D, while the Protein B dimer remains free. This shows that non-phosphorylated Protein A has a higher affinity for C and D than Protein B does.")

    print("\n--- Step 3: Analysis of Experiment 3 (Proteins A, B, C, D + Kinase) ---")
    print("This step investigates the effect of phosphorylation.")
    
    peak1_exp3 = 25
    peak2_exp3 = 40
    peak3_exp3 = 460

    print(f"Observed peaks are at {peak1_exp3} kDa, {peak2_exp3} kDa, and {peak3_exp3} kDa.")
    print(f"The {peak2_exp3} kDa peak is the free kinase.")
    print(f"The {peak1_exp3} kDa peak is crucial: it matches the theoretical mass of a Protein A MONOMER. This strongly implies the Protein A dimer dissociated and that Protein A is the substrate for the kinase.")
    
    complex_BCD_mass = mass_B_native + mass_C_mono + mass_D_mono
    print(f"The {peak3_exp3} kDa peak corresponds to a complex of Protein B dimer + Protein C + Protein D.")
    print(f"Equation: {mass_B_native} + {mass_C_mono} + {mass_D_mono} = {complex_BCD_mass}")
    print("Conclusion: Phosphorylation of Protein A causes it to dissociate into monomers and lose its affinity for C and D. This allows the Protein B dimer to bind to C and D instead.")

    print("\n--- Step 4: Analysis of Experiment 4 (Dephosphorylation) ---")
    print("This step confirms the role of phosphorylation.")

    peak1_exp4 = 50
    peak2_exp4 = 460

    print(f"After removing the kinase and dephosphorylating the sample, the observed peaks are {peak1_exp4} kDa and {peak2_exp4} kDa.")
    print(f"The {peak1_exp4} kDa peak shows that dephosphorylated Protein A monomers re-form a dimer.")
    print(f"Equation: {mass_A_mono} * 2 = {mass_A_mono * 2}")
    print(f"The {peak2_exp4} kDa peak shows the Protein B-C-D complex remains stable and does not get displaced by the now active (non-phosphorylated) Protein A dimer.")
    print("Conclusion: While Protein A has a higher affinity for C and D (Exp 2), the pre-formed B-C-D complex is kinetically stable and does not easily dissociate.")
    
    print("\n--- Step 5: Final Conclusion Evaluation ---")
    print("Based on the analysis, we have the following facts:")
    print("1. Non-phosphorylated Protein A has a higher affinity for C+D than Protein B does (from Exp 2).")
    print("2. Phosphorylation of Protein A reduces its affinity for C+D and causes its dimer to dissociate (from Exp 3).")
    print("3. Protein B always exists as a dimer. Protein A exists as a dimer or a monomer. Proteins C and D are monomers.")
    
    print("\nEvaluating Answer Choice G:")
    print("'G. Protein B never has a higher affinity for protein C or D than nonphosphorylated protein A, protein B always exists as a dimer in opposition to proteins A, C, and D.'")
    print("- 'Protein B never has a higher affinity...than nonphosphorylated protein A': This is supported by Experiment 2, which shows A outcompeting B. The result in Experiment 4 is a kinetic effect, not a reflection of higher affinity.")
    print("- 'protein B always exists as a dimer...': This is true based on Experiment 1.")
    print("- '...in opposition to proteins A, C, and D.': This is also true. B is always a dimer, whereas A can be a monomer or dimer, and C and D are always monomers. B's state is distinct from the others.")
    print("\nTherefore, Answer G is the most correct statement based on the provided data.")

solve_protein_puzzle()
print("<<<G>>>")