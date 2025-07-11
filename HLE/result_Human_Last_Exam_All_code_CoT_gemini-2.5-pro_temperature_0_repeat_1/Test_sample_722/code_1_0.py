def analyze_protein_interactions():
    """
    Analyzes the provided experimental data to determine the correct conclusion.
    """
    # Theoretical masses (in kDa)
    theo_A = 25
    theo_B = 150
    theo_C = 60
    theo_D = 100
    kinase = 40

    print("--- Step-by-step Analysis ---")

    # --- Experiment 1: Individual Proteins ---
    print("\nStep 1: Analysis of Experiment 1 (Individual Proteins)")
    exp1_A = 50
    exp1_B = 300
    exp1_C = 60
    exp1_D = 100
    
    # Determine oligomeric states
    state_A = exp1_A / theo_A
    state_B = exp1_B / theo_B
    
    print(f"Protein A: Experimental mass is {exp1_A} kDa. Since theoretical mass is {theo_A} kDa, it exists as a homodimer ({int(state_A)} * {theo_A} = {exp1_A}).")
    print(f"Protein B: Experimental mass is {exp1_B} kDa. Since theoretical mass is {theo_B} kDa, it exists as a homodimer ({int(state_B)} * {theo_B} = {exp1_B}).")
    print(f"Protein C: Experimental mass ({exp1_C} kDa) matches theoretical mass ({theo_C} kDa), so it's a monomer.")
    print(f"Protein D: Experimental mass ({exp1_D} kDa) matches theoretical mass ({theo_D} kDa), so it's a monomer.")
    
    # Store the native masses for convenience
    native_A_dimer = exp1_A
    native_B_dimer = exp1_B
    native_C = exp1_C
    native_D = exp1_D

    # --- Experiment 2: All Four Proteins Mixed ---
    print("\nStep 2: Analysis of Experiment 2 (A+B+C+D)")
    exp2_peak1 = 300
    exp2_peak2 = 210
    
    print(f"Peak 1 at {exp2_peak1} kDa corresponds to the free Protein B dimer.")
    print(f"Peak 2 is at {exp2_peak2} kDa. Let's check the combination of the other proteins:")
    complex_ACD = native_A_dimer + native_C + native_D
    print(f"Equation: Protein A (dimer) + Protein C + Protein D = {native_A_dimer} + {native_C} + {native_D} = {complex_ACD} kDa.")
    print("Conclusion: This matches Peak 2. Therefore, the A(dimer)-C-D complex forms, implying nonphosphorylated Protein A has a higher affinity for C and D than Protein B does.")

    # --- Experiment 3: Mixed with Kinase ---
    print("\nStep 3: Analysis of Experiment 3 (A+B+C+D + Kinase)")
    exp3_peak1 = 25
    exp3_peak2 = 40
    exp3_peak3 = 460
    
    print(f"Peak 1 at {exp3_peak1} kDa matches the theoretical mass of monomeric Protein A ({theo_A} kDa).")
    print(f"Peak 2 at {exp3_peak2} kDa matches the mass of the free kinase.")
    print("This suggests the kinase phosphorylates Protein A, causing its dimer to dissociate into monomers.")
    print(f"Peak 3 is at {exp3_peak3} kDa. Let's check the combination of the remaining proteins:")
    complex_BCD = native_B_dimer + native_C + native_D
    print(f"Equation: Protein B (dimer) + Protein C + Protein D = {native_B_dimer} + {native_C} + {native_D} = {complex_BCD} kDa.")
    print("Conclusion: This matches Peak 3. Phosphorylated, monomeric Protein A loses its affinity for C and D, allowing the B(dimer)-C-D complex to form.")

    # --- Experiment 4: Dephosphorylated Sample ---
    print("\nStep 4: Analysis of Experiment 4 (Dephosphorylated A, B, C, D)")
    exp4_peak1 = 50
    exp4_peak2 = 460
    
    print(f"Peak 1 at {exp4_peak1} kDa corresponds to the re-formed Protein A dimer.")
    print(f"Peak 2 at {exp4_peak2} kDa corresponds to the B-C-D complex, which remains intact.")
    print(f"Equation: Protein B (dimer) + Protein C + Protein D = {native_B_dimer} + {native_C} + {native_D} = {exp4_peak2} kDa.")
    print("Conclusion: The B-C-D complex is stable and does not dissociate even when the A-dimer is present again.")

    # --- Final Evaluation ---
    print("\nStep 5: Evaluating the Answer Choices based on the analysis")
    print("A, C, H are incorrect because the evidence points to phosphorylation of Protein A, not B.")
    print("B, D are incorrect because phosphorylation of Protein A *decreases* its affinity for C and D.")
    print("E is incorrect because Protein A exists as a monomer in Experiment 3.")
    print("F is incorrect because nonphosphorylated Protein A has a higher affinity for C+D than B does (Exp 2), and A is not always a monomer.")
    print("G is correct because:")
    print("  - Exp 2 shows Protein B has a lower affinity for C+D than nonphosphorylated Protein A.")
    print("  - All experiments show Protein B's mass contribution is 300 kDa, meaning it is always a dimer.")
    print("  - Protein A can be a monomer or dimer, while C and D are monomers, so B is unique in its constant dimeric state.")

analyze_protein_interactions()
<<<G>>>