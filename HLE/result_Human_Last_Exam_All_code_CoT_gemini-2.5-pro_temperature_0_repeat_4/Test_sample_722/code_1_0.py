def solve_protein_puzzle():
    """
    Analyzes the SEC-MALS data to determine the nature of protein interactions.
    """
    # Theoretical masses (kDa)
    mass_A_mono = 25
    mass_B_mono = 150
    mass_C_mono = 60
    mass_D_mono = 100
    mass_kinase = 40

    print("Step 1: Analysis of Experiment 1 (Individual Proteins)")
    # Experimental masses from Exp 1
    exp1_A = 50
    exp1_B = 300
    exp1_C = 60
    exp1_D = 100
    
    # Determine oligomeric states
    state_A = exp1_A / mass_A_mono
    state_B = exp1_B / mass_B_mono
    state_C = exp1_C / mass_C_mono
    state_D = exp1_D / mass_D_mono

    print(f"Protein A: Experimental mass = {exp1_A} kDa. Theoretical monomer = {mass_A_mono} kDa. Oligomeric state = {int(state_A)} (Dimer).")
    print(f"Protein B: Experimental mass = {exp1_B} kDa. Theoretical monomer = {mass_B_mono} kDa. Oligomeric state = {int(state_B)} (Dimer).")
    print(f"Protein C: Experimental mass = {exp1_C} kDa. Theoretical monomer = {mass_C_mono} kDa. Oligomeric state = {int(state_C)} (Monomer).")
    print(f"Protein D: Experimental mass = {exp1_D} kDa. Theoretical monomer = {mass_D_mono} kDa. Oligomeric state = {int(state_D)} (Monomer).")
    print("-" * 30)

    # Define the masses of the stable forms of the proteins
    mass_A_dimer = exp1_A
    mass_B_dimer = exp1_B

    print("Step 2: Analysis of Experiment 2 (Mixture of A, B, C, D)")
    # Experimental peaks from Exp 2
    exp2_peak1 = 300
    exp2_peak2 = 210
    
    # Interpret peaks
    print(f"Observed Peak 1: {exp2_peak1} kDa. This corresponds to the mass of the Protein B dimer ({mass_B_dimer} kDa), which is unbound.")
    
    complex_A2CD = mass_A_dimer + mass_C_mono + mass_D_mono
    print(f"Observed Peak 2: {exp2_peak2} kDa. Let's test combinations:")
    print(f"Protein A Dimer ({mass_A_dimer}) + Protein C ({mass_C_mono}) + Protein D ({mass_D_mono}) = {complex_A2CD} kDa.")
    print("This matches Peak 2. So, a complex of A2CD is formed.")
    print("Conclusion: In the absence of phosphorylation, the A2CD complex is favored, and Protein B remains free.")
    print("-" * 30)

    print("Step 3: Analysis of Experiment 3 (Mixture + Kinase)")
    # Experimental peaks from Exp 3
    exp3_peak1 = 25
    exp3_peak2 = 40
    exp3_peak3 = 460

    print(f"Observed Peak 1: {exp3_peak1} kDa. This matches the theoretical mass of monomeric Protein A ({mass_A_mono} kDa).")
    print(f"Observed Peak 2: {exp3_peak2} kDa. This matches the mass of the kinase ({mass_kinase} kDa).")
    
    complex_B2CD = mass_B_dimer + mass_C_mono + mass_D_mono
    print(f"Observed Peak 3: {exp3_peak3} kDa. Let's test combinations:")
    print(f"Protein B Dimer ({mass_B_dimer}) + Protein C ({mass_C_mono}) + Protein D ({mass_D_mono}) = {complex_B2CD} kDa.")
    print("This matches Peak 3. So, a complex of B2CD is formed.")
    print("Conclusion: The presence of the kinase causes the A2CD complex to dissociate and a new B2CD complex to form. Protein A dissociates into a monomer.")
    print("-" * 30)
    
    print("Step 4: Analysis of Experiment 4 (Dephosphorylation)")
    print("The sample from Exp 3 was treated with phosphatase, which removed the phosphate from Protein A.")
    # Experimental peaks from Exp 4
    exp4_peak1 = 50
    exp4_peak2 = 460
    print(f"Observed Peak 1: {exp4_peak1} kDa. This is the mass of the Protein A dimer ({mass_A_dimer} kDa).")
    print(f"Observed Peak 2: {exp4_peak2} kDa. This is the mass of the B2CD complex ({complex_B2CD} kDa).")
    print("Conclusion: Dephosphorylation allows Protein A to re-form a dimer, but the B2CD complex is stable and does not dissociate.")
    print("-" * 30)

    print("Step 5: Evaluating the Answer Choices")
    print("Hypothesis: To explain the results of Exp 3, the kinase must have caused the observed switch. The dissociation of Protein A into a monomer suggests it was phosphorylated, decreasing its affinity for itself and for C and D. The formation of the B2CD complex suggests phosphorylation increased Protein B's affinity for C and D. The most complete hypothesis is that the kinase acts on both A and B.")
    print("\nLet's evaluate Choice C based on this hypothesis:")
    print("Choice C: 'Phosphorylation of protein B increases its affinity for protein C or D and dephosphorylated protein A always has a higher affinity for protein C or D than protein B. Protein A exists as a monomer or homodimer.'")
    
    print("\nClause 1: 'Phosphorylation of protein B increases its affinity for protein C or D'")
    print("This is consistent with the formation of the B2CD complex in Exp 3 only after the kinase is added. So, this clause is TRUE.")
    
    print("\nClause 2: 'dephosphorylated protein A always has a higher affinity for protein C or D than protein B'")
    print("This is consistent with Exp 2, where the A2CD complex forms preferentially over a B2CD complex. So, this clause is TRUE.")

    print("\nClause 3: 'Protein A exists as a monomer or homodimer.'")
    print("This is consistent with Exp 1/4 (dimer, 50 kDa) and Exp 3 (monomer, 25 kDa). So, this clause is TRUE.")

    print("\nFinal Conclusion: All three clauses of choice C are consistent with the experimental data under the most complete hypothesis.")
    print("The key interaction equations are:")
    print(f"Equation (Dimerization A): {mass_A_mono} + {mass_A_mono} -> {mass_A_dimer}")
    print(f"Equation (Dimerization B): {mass_B_mono} + {mass_B_mono} -> {mass_B_dimer}")
    print(f"Equation (Baseline Complex): {mass_A_dimer} + {mass_C_mono} + {mass_D_mono} -> {complex_A2CD}")
    print(f"Equation (Phosphorylated Complex): {mass_B_dimer} + {mass_C_mono} + {mass_D_mono} -> {complex_B2CD}")

solve_protein_puzzle()
print("<<<C>>>")