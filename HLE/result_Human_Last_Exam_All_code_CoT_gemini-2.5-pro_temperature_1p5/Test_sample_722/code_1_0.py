def analyze_protein_interactions():
    """
    Analyzes SEC-MALS experimental data to deduce protein interactions
    and evaluate given conclusions.
    """
    # Theoretical and observed masses (in kDa)
    proteins = {
        'A': {'theoretical': 25, 'exp1': 50},
        'B': {'theoretical': 150, 'exp1': 300},
        'C': {'theoretical': 60, 'exp1': 60},
        'D': {'theoretical': 100, 'exp1': 100},
        'Kinase': {'mass': 40}
    }

    print("--- Step 1: Analysis of Experiment 1 (Oligomeric States) ---\n")
    
    # Determine oligomeric states
    protein_a_state = proteins['A']['exp1'] / proteins['A']['theoretical']
    protein_b_state = proteins['B']['exp1'] / proteins['B']['theoretical']
    
    print(f"Protein A: Theoretical={proteins['A']['theoretical']}kDa, Observed={proteins['A']['exp1']}kDa -> Forms a Dimer ({int(protein_a_state)} * {proteins['A']['theoretical']} = {proteins['A']['exp1']})")
    print(f"Protein B: Theoretical={proteins['B']['theoretical']}kDa, Observed={proteins['B']['exp1']}kDa -> Forms a Dimer ({int(protein_b_state)} * {proteins['B']['theoretical']} = {proteins['B']['exp1']})")
    print(f"Protein C: Theoretical={proteins['C']['theoretical']}kDa, Observed={proteins['C']['exp1']}kDa -> Exists as a Monomer")
    print(f"Protein D: Theoretical={proteins['D']['theoretical']}kDa, Observed={proteins['D']['exp1']}kDa -> Exists as a Monomer\n")

    # Masses of the units that will be used in subsequent experiments
    a_dimer = proteins['A']['exp1']
    b_dimer = proteins['B']['exp1']
    c_monomer = proteins['C']['exp1']
    d_monomer = proteins['D']['exp1']

    print("--- Step 2: Analysis of Experiment 2 (Protein Mixture) ---\n")
    exp2_peak1 = 300
    exp2_peak2 = 210
    
    print(f"Observed peaks: {exp2_peak1} kDa and {exp2_peak2} kDa.")
    print(f"The {exp2_peak1} kDa peak corresponds to the free Protein B dimer ({b_dimer} kDa).")
    complex_2_calc = a_dimer + c_monomer + d_monomer
    print(f"The {exp2_peak2} kDa peak is a complex of: Protein A dimer ({a_dimer} kDa) + Protein C ({c_monomer} kDa) + Protein D ({d_monomer} kDa) = {complex_2_calc} kDa.")
    print("Conclusion: Dephosphorylated Protein A dimer has a higher affinity for Proteins C and D than the Protein B dimer does.\n")
    
    print("--- Step 3: Analysis of Experiment 3 (Mixture with Kinase) ---\n")
    exp3_peak1 = 25
    exp3_peak2 = 40
    exp3_peak3 = 460
    
    print(f"Observed peaks: {exp3_peak1} kDa, {exp3_peak2} kDa, and {exp3_peak3} kDa.")
    print(f"The {exp3_peak2} kDa peak is the free Kinase ({proteins['Kinase']['mass']} kDa).")
    print(f"The {exp3_peak1} kDa peak is a Protein A monomer ({proteins['A']['theoretical']} kDa). This indicates phosphorylation by the kinase caused the Protein A dimer to dissociate.")
    complex_3_calc = b_dimer + c_monomer + d_monomer
    print(f"The {exp3_peak3} kDa peak is a complex of: Protein B dimer ({b_dimer} kDa) + Protein C ({c_monomer} kDa) + Protein D ({d_monomer} kDa) = {complex_3_calc} kDa.")
    print("Conclusion: Phosphorylation of Protein A decreases its affinity for Proteins C and D, allowing Protein B to form a complex with them.\n")
    
    print("--- Step 4: Analysis of Experiment 4 (Dephosphorylation) ---\n")
    exp4_peak1 = 50
    exp4_peak2 = 460
    
    print(f"Observed peaks: {exp4_peak1} kDa and {exp4_peak2} kDa.")
    print(f"The {exp4_peak1} kDa peak is the free Protein A dimer ({a_dimer} kDa), which re-formed after dephosphorylation.")
    print(f"The {exp4_peak2} kDa peak is the stable complex of Protein B, C, and D ({b_dimer} kDa + {c_monomer} kDa + {d_monomer} kDa = {complex_3_calc} kDa) that formed in Experiment 3.")
    print("Conclusion: The Protein B/C/D complex, once formed, is stable and is not displaced by the dephosphorylated Protein A dimer.\n")

    print("--- Step 5: Final Evaluation of Answer Choices ---\n")
    print("A: Incorrect. Protein A is phosphorylated, not B.")
    print("B: Incorrect. Phosphorylation of A *decreases* its affinity for C and D.")
    print("C: Incorrect. Protein A is phosphorylated, not B.")
    print("D: Incorrect. Phosphorylation of A *decreases* its affinity.")
    print("E: Incorrect. Protein A exists as a monomer when phosphorylated.")
    print("F: Incorrect. Protein A exists as a dimer in its native state.")
    print("G: Incorrect. Contradicted by Exp 3/4 and makes false claims on oligomeric states.")
    print("H: Incorrect. Protein A is phosphorylated, not B.")
    print("I: Incorrect as C and H are incorrect.")
    print("\nBased on the analysis, none of the answer choices A through I are fully correct.")

analyze_protein_interactions()
print("<<<J>>>")