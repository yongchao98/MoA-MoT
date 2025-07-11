def analyze_protein_interactions():
    """
    Analyzes the SEC-MALS experimental data to find the most correct conclusion.
    """
    # Define protein masses based on the problem description
    theoretical_mass = {'A': 25, 'B': 150, 'C': 60, 'D': 100}
    kinase_mass = 40

    # --- Introduction ---
    print("Analyzing the SEC-MALS experiments step-by-step to determine the correct conclusion.\n")

    # --- Experiment 1 Analysis ---
    print("--- Analysis of Experiment 1: Individual Proteins ---")
    exp1_masses = {'A': 50, 'B': 300, 'C': 60, 'D': 100}
    oligomeric_state_A = exp1_masses['A'] / theoretical_mass['A']
    oligomeric_state_B = exp1_masses['B'] / theoretical_mass['B']

    print(f"Protein A: Theoretical mass = {theoretical_mass['A']} kDa, Experimental mass = {exp1_masses['A']} kDa. This is a homodimer ({int(oligomeric_state_A)} * {theoretical_mass['A']} = {exp1_masses['A']}).")
    print(f"Protein B: Theoretical mass = {theoretical_mass['B']} kDa, Experimental mass = {exp1_masses['B']} kDa. This is a homodimer ({int(oligomeric_state_B)} * {theoretical_mass['B']} = {exp1_masses['B']}).")
    print(f"Protein C: Theoretical mass = {theoretical_mass['C']} kDa, Experimental mass = {exp1_masses['C']} kDa. This is a monomer.")
    print(f"Protein D: Theoretical mass = {theoretical_mass['D']} kDa, Experimental mass = {exp1_masses['D']} kDa. This is a monomer.")
    print("Conclusion from Exp 1: In their native state, A and B are homodimers, C and D are monomers.\n")

    # Store native state masses for later use
    native_mass = {
        'A_dimer': exp1_masses['A'],
        'B_dimer': exp1_masses['B'],
        'C_monomer': exp1_masses['C'],
        'D_monomer': exp1_masses['D']
    }

    # --- Experiment 2 Analysis ---
    print("--- Analysis of Experiment 2: Mixture of A, B, C, D ---")
    exp2_peaks = [300, 210]
    print(f"Observed peaks at {exp2_peaks[0]} kDa and {exp2_peaks[1]} kDa.")

    # Identify the complexes
    complex_ACD_mass = native_mass['A_dimer'] + native_mass['C_monomer'] + native_mass['D_monomer']
    print(f"The 210 kDa peak corresponds to a complex of Protein A dimer, Protein C, and Protein D.")
    print(f"Calculation: {native_mass['A_dimer']} kDa (A dimer) + {native_mass['C_monomer']} kDa (C) + {native_mass['D_monomer']} kDa (D) = {complex_ACD_mass} kDa.")
    print(f"The 300 kDa peak corresponds to unbound Protein B dimer.")
    print("Conclusion from Exp 2: In a direct competition, nonphosphorylated Protein A has a higher affinity for C and D than Protein B does.\n")


    # --- Experiment 3 Analysis ---
    print("--- Analysis of Experiment 3: Mixture + Kinase ---")
    exp3_peaks = [25, 40, 460]
    print(f"Observed peaks at {exp3_peaks[0]} kDa, {exp3_peaks[1]} kDa, and {exp3_peaks[2]} kDa.")

    print(f"The {exp3_peaks[1]} kDa peak is the unbound kinase.")
    print(f"The {exp3_peaks[0]} kDa peak is the Protein A monomer, indicating its dimer dissociated upon phosphorylation.")

    complex_BCD_mass = native_mass['B_dimer'] + native_mass['C_monomer'] + native_mass['D_monomer']
    print(f"The {exp3_peaks[2]} kDa peak corresponds to a complex of Protein B dimer, Protein C, and Protein D.")
    print(f"Calculation: {native_mass['B_dimer']} kDa (B dimer) + {native_mass['C_monomer']} kDa (C) + {native_mass['D_monomer']} kDa (D) = {complex_BCD_mass} kDa.")
    print("Conclusion from Exp 3: Phosphorylation of Protein A prevents it from binding C and D, allowing Protein B to form a complex.\n")

    # --- Experiment 4 Analysis ---
    print("--- Analysis of Experiment 4: Dephosphorylation of Protein A ---")
    print(f"The prompt confirms Protein A was the phosphorylated protein.")
    print(f"Observed peaks are {exp4_peaks[0]} kDa and {exp4_peaks[1]} kDa.")
    print(f"The {exp4_peaks[0]} kDa peak is the reformed Protein A dimer ({native_mass['A_dimer']} kDa).")
    print(f"The {exp4_peaks[1]} kDa peak is the stable Protein B-C-D complex ({complex_BCD_mass} kDa).")
    print("Conclusion from Exp 4: The B-C-D complex, once formed, is kinetically stable and not displaced by dephosphorylated Protein A.\n")

    # --- Final Evaluation of Answer Choices ---
    print("--- Evaluating the Answer Choices ---")
    print("Choice G states: 'Protein B never has a higher affinity for protein C or D than nonphosphorylated protein A, protein B always exists as a dimer in opposition to proteins A, C, and D.'")
    print(" - Part 1 ('Protein B never has a higher affinity...'): This is supported by Experiment 2, where nonphosphorylated A outcompetes B.")
    print(" - Part 2 ('protein B always exists as a dimer...'): This is supported by all experiments. B is always a dimer, while A's state changes and C/D are monomers.")
    print("\nTherefore, Choice G is the most accurate conclusion.")

analyze_protein_interactions()