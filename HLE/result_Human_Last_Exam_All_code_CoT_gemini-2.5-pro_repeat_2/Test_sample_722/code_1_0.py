def analyze_protein_interactions():
    """
    Analyzes the SEC-MALS experimental data to determine the correct conclusion.
    """
    # Theoretical masses of the monomers
    theo_mass = {'A': 25, 'B': 150, 'C': 60, 'D': 100}
    kinase_mass = 40

    print("--- Analysis of Experiment 1: Determining Oligomeric States ---")
    exp1_mass = {'A': 50, 'B': 300, 'C': 60, 'D': 100}
    
    # Calculate and print oligomeric states
    a_dimer_mass = exp1_mass['A']
    b_dimer_mass = exp1_mass['B']
    c_mono_mass = exp1_mass['C']
    d_mono_mass = exp1_mass['D']

    print(f"Protein A: Observed {a_dimer_mass} kDa / Theoretical {theo_mass['A']} kDa = {int(a_dimer_mass / theo_mass['A'])}. It's a dimer.")
    print(f"Protein B: Observed {b_dimer_mass} kDa / Theoretical {theo_mass['B']} kDa = {int(b_dimer_mass / theo_mass['B'])}. It's a dimer.")
    print(f"Protein C: Observed {c_mono_mass} kDa / Theoretical {theo_mass['C']} kDa = {int(c_mono_mass / theo_mass['C'])}. It's a monomer.")
    print(f"Protein D: Observed {d_mono_mass} kDa / Theoretical {theo_mass['D']} kDa = {int(d_mono_mass / theo_mass['D'])}. It's a monomer.")

    print("\n--- Analysis of Experiment 2: Baseline Interactions (No Kinase) ---")
    exp2_peak1 = 300
    exp2_peak2 = 210
    
    a_c_d_complex = a_dimer_mass + c_mono_mass + d_mono_mass
    print(f"Peak 1 at {exp2_peak1} kDa matches the free Protein B dimer.")
    print(f"Peak 2 at {exp2_peak2} kDa matches the sum of Protein A dimer ({a_dimer_mass}) + C ({c_mono_mass}) + D ({d_mono_mass}), which is {a_c_d_complex} kDa.")
    print("Conclusion: Non-phosphorylated Protein A outcompetes Protein B for binding to C and D.")

    print("\n--- Analysis of Experiment 3: Effect of Phosphorylation ---")
    exp3_peak1 = 25
    exp3_peak2 = 40
    exp3_peak3 = 460
    
    b_c_d_complex = b_dimer_mass + c_mono_mass + d_mono_mass
    print(f"Peak 1 at {exp3_peak1} kDa matches the Protein A monomer ({theo_mass['A']} kDa). The A-C-D complex and A-A dimer broke apart.")
    print(f"Peak 2 at {exp3_peak2} kDa matches the free kinase.")
    print(f"Peak 3 at {exp3_peak3} kDa matches the sum of Protein B dimer ({b_dimer_mass}) + C ({c_mono_mass}) + D ({d_mono_mass}), which is {b_c_d_complex} kDa.")
    print("Conclusion: Phosphorylation of Protein A disrupts its interactions, allowing B-C-D to form.")

    print("\n--- Analysis of Experiment 4: Effect of Dephosphorylation ---")
    exp4_peak1 = 50
    exp4_peak2 = 460
    print(f"Peak 1 at {exp4_peak1} kDa shows Protein A has re-formed its dimer ({a_dimer_mass} kDa).")
    print(f"Peak 2 at {exp4_peak2} kDa shows the B-C-D complex ({b_c_d_complex} kDa) remains intact.")
    print("Conclusion: The B-C-D complex is highly stable once formed.")
    
    print("\n--- Final Evaluation of Answer Choice G ---")
    print("Statement G: 'Protein B never has a higher affinity for protein C or D than nonphosphorylated protein A, protein B always exists as a dimer in opposition to proteins A, C, and D.'")
    print("1. 'Protein B never has a higher affinity... than nonphosphorylated protein A': Correct. Exp 2 shows A-C-D forms, leaving B free.")
    print("2. 'protein B always exists as a dimer...': Correct. Exp 1 shows B is a dimer, while A can be a monomer/dimer and C/D are monomers.")
    print("Both clauses of statement G are fully supported by the data.")

analyze_protein_interactions()