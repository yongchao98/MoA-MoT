def solve_protein_interaction_puzzle():
    """
    Analyzes SEC-MALS data to determine protein interactions and select the correct conclusion.
    """

    # --- Step 1: Define Theoretical and Experimental Masses ---
    print("--- Analysis Start ---")
    print("\nStep 1: Define Protein Masses and Analyze Experiment 1 (Individual Proteins)")

    theoretical_masses = {'A': 25, 'B': 150, 'C': 60, 'D': 100, 'Kinase': 40}
    exp1_masses = {'A': 50, 'B': 300, 'C': 60, 'D': 100}

    # Calculate and print oligomeric states
    oligomeric_states = {}
    for protein, exp_mass in exp1_masses.items():
        theo_mass = theoretical_masses[protein]
        state = round(exp_mass / theo_mass)
        state_str = "monomer" if state == 1 else "dimer" if state == 2 else f"{state}-mer"
        oligomeric_states[protein] = {'mass': exp_mass, 'state': state_str}
        print(f"Protein {protein}: Theoretical mass = {theo_mass} kDa, Experimental mass = {exp_mass} kDa. "
              f"Calculation: {exp_mass} / {theo_mass} = {state}. Conclusion: Protein {protein} is a {state_str}.")

    # --- Step 2: Analyze Experiment 2 (Mixture of A, B, C, D) ---
    print("\nStep 2: Analyze Experiment 2 (Mixture of A, B, C, D)")
    exp2_peaks = [300, 210]
    print(f"Observed peaks at {exp2_peaks[0]} kDa and {exp2_peaks[1]} kDa.")

    a_mass = oligomeric_states['A']['mass']
    b_mass = oligomeric_states['B']['mass']
    c_mass = oligomeric_states['C']['mass']
    d_mass = oligomeric_states['D']['mass']

    print(f"The {exp2_peaks[0]} kDa peak corresponds to the mass of the Protein B dimer ({b_mass} kDa). This means Protein B is not in a complex.")
    
    complex_acd_mass = a_mass + c_mass + d_mass
    print(f"Let's test if the remaining proteins (A, C, D) form the second peak:")
    print(f"Mass of A ({oligomeric_states['A']['state']}) + C ({oligomeric_states['C']['state']}) + D ({oligomeric_states['D']['state']}) = {a_mass} + {c_mass} + {d_mass} = {complex_acd_mass} kDa.")
    print(f"This matches the {exp2_peaks[1]} kDa peak.")
    print("Conclusion: Non-phosphorylated Protein A (dimer) forms a complex with C and D. Protein B is excluded. This indicates Affinity(A for C+D) > Affinity(B for C+D).")

    # --- Step 3: Analyze Experiment 3 (Mixture + Kinase) ---
    print("\nStep 3: Analyze Experiment 3 (Mixture + Kinase)")
    exp3_peaks = [25, 40, 460]
    print(f"Observed peaks at {exp3_peaks[0]} kDa, {exp3_peaks[1]} kDa, and {exp3_peaks[2]} kDa.")

    print(f"The {exp3_peaks[0]} kDa peak corresponds to the theoretical mass of monomeric Protein A ({theoretical_masses['A']} kDa). This suggests phosphorylation by the kinase caused the Protein A dimer to dissociate into monomers.")
    print(f"The {exp3_peaks[1]} kDa peak corresponds to the mass of the free kinase ({theoretical_masses['Kinase']} kDa).")

    complex_bcd_mass = b_mass + c_mass + d_mass
    print(f"Let's test the composition of the large complex:")
    print(f"Mass of B ({oligomeric_states['B']['state']}) + C ({oligomeric_states['C']['state']}) + D ({oligomeric_states['D']['state']}) = {b_mass} + {c_mass} + {d_mass} = {complex_bcd_mass} kDa.")
    print(f"This matches the {exp3_peaks[2]} kDa peak.")
    print("Conclusion: Phosphorylation of Protein A causes it to dissociate and leave the complex. Protein B then binds to C and D.")

    # --- Step 4: Analyze Experiment 4 (Dephosphorylation of A) ---
    print("\nStep 4: Analyze Experiment 4 (Dephosphorylation of A)")
    exp4_peaks = [50, 460]
    print(f"Observed peaks at {exp4_peaks[0]} kDa and {exp4_peaks[1]} kDa.")

    print(f"The {exp4_peaks[0]} kDa peak corresponds to the mass of the Protein A dimer ({a_mass} kDa). This shows dephosphorylation allowed Protein A to re-form its dimer.")
    print(f"The {exp4_peaks[1]} kDa peak is the same B-C-D complex from Experiment 3 ({complex_bcd_mass} kDa), which remains intact.")
    print("Conclusion: Even when Protein A is dephosphorylated and dimerizes, it does not displace Protein B from the stable B-C-D complex once it has formed.")
    
    # --- Step 5: Evaluate Answer Choices and Final Conclusion ---
    print("\nStep 5: Final Conclusion based on Analysis")
    print("Based on the experiments:")
    print("1. From Exp 2, non-phosphorylated Protein A has a higher affinity for C+D than Protein B does, as A forms the complex and B is left free.")
    print("2. From all experiments, Protein B is consistently a dimer (300 kDa).")
    print("3. Protein A can be a dimer (Exp 1, 2, 4) or a monomer (Exp 3). Proteins C and D are always monomers.")
    print("Therefore, Protein B is the only protein that is *always* a dimer.")
    
    print("\nEvaluating choice G: 'Protein B never has a higher affinity for protein C or D than nonphosphorylated protein A, protein B always exists as a dimer in opposition to proteins A, C, and D.'")
    print("- 'Protein B never has a higher affinity... than nonphosphorylated protein A' -> Correct, as shown in Exp 2.")
    print("- 'protein B always exists as a dimer' -> Correct, its mass is always 300 kDa.")
    print("- 'in opposition to proteins A, C, and D' -> Correct, A can be a monomer, C and D are monomers.")
    print("Statement G is fully supported by the data.")

# Execute the analysis
solve_protein_interaction_puzzle()

# Final Answer
print("\n<<<G>>>")