def analyze_sec_mals_data():
    """
    Analyzes the provided SEC-MALS experimental data to determine the correct conclusion.
    The analysis is printed step-by-step.
    """

    # --- Protein masses (kDa) ---
    theoretical_mass = {'A': 25, 'B': 150, 'C': 60, 'D': 100}
    exp1_mass = {'A': 50, 'B': 300, 'C': 60, 'D': 100}

    print("Step 1: Determine the native state of each protein from Experiment 1.")
    print(f"- Theoretical masses (kDa): A={theoretical_mass['A']}, B={theoretical_mass['B']}, C={theoretical_mass['C']}, D={theoretical_mass['D']}")
    print(f"- Experimental masses from Exp 1 (kDa): A={exp1_mass['A']}, B={exp1_mass['B']}, C={exp1_mass['C']}, D={exp1_mass['D']}\n")

    # Calculate and print oligomeric states
    a_state = exp1_mass['A'] / theoretical_mass['A']
    b_state = exp1_mass['B'] / theoretical_mass['B']
    c_state = exp1_mass['C'] / theoretical_mass['C']
    d_state = exp1_mass['D'] / theoretical_mass['D']

    print(f"Oligomeric state calculations:")
    print(f"- Protein A: {exp1_mass['A']} kDa / {theoretical_mass['A']} kDa = {a_state:.1f}. This indicates Protein A is a homodimer (A2).")
    print(f"- Protein B: {exp1_mass['B']} kDa / {theoretical_mass['B']} kDa = {b_state:.1f}. This indicates Protein B is a homodimer (B2).")
    print(f"- Protein C: {exp1_mass['C']} kDa / {theoretical_mass['C']} kDa = {c_state:.1f}. This indicates Protein C is a monomer (C).")
    print(f"- Protein D: {exp1_mass['D']} kDa / {theoretical_mass['D']} kDa = {d_state:.1f}. This indicates Protein D is a monomer (D).\n")

    # --- Native state masses for complex calculation ---
    a_dimer_mass = exp1_mass['A']
    b_dimer_mass = exp1_mass['B']
    c_mono_mass = exp1_mass['C']
    d_mono_mass = exp1_mass['D']

    print("Step 2: Analyze the complexes formed in Experiment 2 (A, B, C, D mixed).")
    exp2_peaks = [300, 210]
    print(f"- Peaks observed at {exp2_peaks[0]} kDa and {exp2_peaks[1]} kDa.")
    print(f"- The {exp2_peaks[0]} kDa peak corresponds to the free Protein B dimer (B2).")
    print(f"- To identify the {exp2_peaks[1]} kDa peak, let's test a combination of the other proteins: Protein A dimer, Protein C monomer, and Protein D monomer.")
    complex_acd_mass = a_dimer_mass + c_mono_mass + d_mono_mass
    print(f"- Mass of A2 + C + D = {a_dimer_mass} + {c_mono_mass} + {d_mono_mass} = {complex_acd_mass} kDa.")
    print("- This perfectly matches the 210 kDa peak, forming an (A2)CD complex.")
    print("Conclusion from Exp 2: In a competitive environment, the nonphosphorylated A2 dimer binds to C and D, while the B2 dimer remains free. This indicates that A2 has a higher affinity for the CD complex than B2 does.\n")


    print("Step 3: Analyze the complexes formed in Experiment 3 (A, B, C, D + Kinase).")
    exp3_peaks = [25, 40, 460]
    print(f"- Peaks observed at {exp3_peaks[0]} kDa, {exp3_peaks[1]} kDa, and {exp3_peaks[2]} kDa.")
    print(f"- The {exp3_peaks[0]} kDa peak is the mass of monomeric Protein A ({theoretical_mass['A']} kDa). The presence of the kinase suggests this is phosphorylated A (pA) that has dissociated from its dimer.")
    print(f"- The {exp3_peaks[1]} kDa peak corresponds to the free kinase.")
    print(f"- To identify the {exp3_peaks[2]} kDa peak, let's test a combination of the remaining proteins: Protein B dimer, Protein C monomer, and Protein D monomer.")
    complex_bcd_mass = b_dimer_mass + c_mono_mass + d_mono_mass
    print(f"- Mass of B2 + C + D = {b_dimer_mass} + {c_mono_mass} + {d_mono_mass} = {complex_bcd_mass} kDa.")
    print("- This perfectly matches the 460 kDa peak, forming a (B2)CD complex.")
    print("Conclusion from Exp 3: Phosphorylation of Protein A causes it to dissociate into a monomer and release the CD complex. The free CD complex is then bound by the B2 dimer. This means phosphorylation DECREASES Protein A's affinity for C and D.\n")
    
    print("Step 4: Evaluate Answer Choice G based on the analysis.")
    print("Answer G states: 'Protein B never has a higher affinity for protein C or D than nonphosphorylated protein A, protein B always exists as a dimer in opposition to proteins A, C, and D.'\n")
    print("- Evaluation of Part 1: 'Protein B never has a higher affinity for protein C or D than nonphosphorylated protein A.'")
    print("  - From Step 2, we concluded that when competing, nonphosphorylated A2 binds to CD, indicating a higher affinity than B2. This part of the statement is supported by the data.\n")
    print("- Evaluation of Part 2: 'protein B always exists as a dimer in opposition to proteins A, C, and D.'")
    print("  - Protein B is observed as a B2 dimer (300 kDa) in Exp 1 & 2, and as part of the (B2)CD complex in Exp 3 & 4. The B2 dimer unit is always present.")
    print("  - Protein A is observed as a dimer (A2, 50 kDa) in Exp 1 & 4, and as a monomer (pA, 25 kDa) in Exp 3.")
    print("  - Proteins C and D consistently act as monomers.")
    print("  - Therefore, B is the only protein that 'always exists as a dimer' unit throughout the experiments, distinguishing it from A, C, and D. This part of the statement is also supported by the data.\n")

    print("Final Conclusion: Both parts of answer choice G are fully supported by the experimental results.")


# Execute the analysis function
analyze_sec_mals_data()
print("<<<G>>>")
