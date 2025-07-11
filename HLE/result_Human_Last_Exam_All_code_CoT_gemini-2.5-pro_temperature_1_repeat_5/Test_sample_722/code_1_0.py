def analyze_protein_interactions():
    """
    Analyzes SEC-MALS data to determine protein interactions and evaluates conclusions.
    """
    # Theoretical masses based on amino acid sequence
    protein_masses = {'A': 25, 'B': 150, 'C': 60, 'D': 100, 'Kinase': 40}

    print("--- Analysis of Protein-Protein Interaction Experiments ---")
    print(f"Theoretical masses (kDa): A={protein_masses['A']}, B={protein_masses['B']}, C={protein_masses['C']}, D={protein_masses['D']}")
    print("-" * 60)

    # --- Step 1: Analysis of Experiment 1 (Individual Proteins) ---
    print("Step 1: Analyzing Experiment 1 (Individual Proteins)")
    exp1_obs = {'A': 50, 'B': 300, 'C': 60, 'D': 100}
    
    oligomer_A = exp1_obs['A'] / protein_masses['A']
    oligomer_B = exp1_obs['B'] / protein_masses['B']
    oligomer_C = exp1_obs['C'] / protein_masses['C']
    oligomer_D = exp1_obs['D'] / protein_masses['D']

    print(f"Protein A: Observed mass {exp1_obs['A']} kDa / Theoretical {protein_masses['A']} kDa = {oligomer_A:.0f}. Conclusion: Protein A is a homodimer.")
    print(f"Protein B: Observed mass {exp1_obs['B']} kDa / Theoretical {protein_masses['B']} kDa = {oligomer_B:.0f}. Conclusion: Protein B is a homodimer.")
    print(f"Protein C: Observed mass {exp1_obs['C']} kDa / Theoretical {protein_masses['C']} kDa = {oligomer_C:.0f}. Conclusion: Protein C is a monomer.")
    print(f"Protein D: Observed mass {exp1_obs['D']} kDa / Theoretical {protein_masses['D']} kDa = {oligomer_D:.0f}. Conclusion: Protein D is a monomer.")
    
    # Store native state masses for later use
    native_masses = {'A': 50, 'B': 300, 'C': 60, 'D': 100}
    print(f"\nSummary of native states: A=homodimer ({native_masses['A']} kDa), B=homodimer ({native_masses['B']} kDa), C=monomer ({native_masses['C']} kDa), D=monomer ({native_masses['D']} kDa)")
    print("-" * 60)

    # --- Step 2: Analysis of Experiment 2 (Mixture of A, B, C, D) ---
    print("Step 2: Analyzing Experiment 2 (Mixture of A, B, C, D)")
    exp2_peaks = [300, 210]
    print(f"Observed peaks at {exp2_peaks[0]} kDa and {exp2_peaks[1]} kDa.")
    print(f"The {exp2_peaks[0]} kDa peak corresponds to the Protein B homodimer ({native_masses['B']} kDa), which is unbound.")
    print(f"The {exp2_peaks[1]} kDa peak is a complex of the remaining proteins: A (dimer), C (monomer), and D (monomer).")
    print(f"Equation: {native_masses['A']} (A_dimer) + {native_masses['C']} (C_monomer) + {native_masses['D']} (D_monomer) = {native_masses['A'] + native_masses['C'] + native_masses['D']} kDa.")
    print("Conclusion for Exp 2: The A(dimer)-C-D complex forms. This means the affinity of non-phosphorylated Protein A for C+D is HIGHER than that of Protein B.")
    print("-" * 60)

    # --- Step 3: Analysis of Experiment 3 (Mixture + Kinase) ---
    print("Step 3: Analyzing Experiment 3 (Mixture + Kinase)")
    exp3_peaks = [25, 40, 460]
    print(f"Observed peaks at {exp3_peaks[0]} kDa, {exp3_peaks[1]} kDa, and {exp3_peaks[2]} kDa.")
    print(f"The {exp3_peaks[1]} kDa peak is the free Kinase ({protein_masses['Kinase']} kDa).")
    print(f"The {exp3_peaks[0]} kDa peak is the mass of monomeric Protein A ({protein_masses['A']} kDa). The kinase must have phosphorylated Protein A, causing its dimer to dissociate.")
    print(f"The {exp3_peaks[2]} kDa peak is a new complex of the remaining proteins: B (dimer), C (monomer), and D (monomer).")
    print(f"Equation: {native_masses['B']} (B_dimer) + {native_masses['C']} (C_monomer) + {native_masses['D']} (D_monomer) = {native_masses['B'] + native_masses['C'] + native_masses['D']} kDa.")
    print("Conclusion for Exp 3: Phosphorylation of Protein A DECREASES its affinity for C+D, allowing the B-C-D complex to form instead.")
    print("-" * 60)

    # --- Step 4: Analysis of Experiment 4 (Dephosphorylation of A) ---
    print("Step 4: Analyzing Experiment 4 (Dephosphorylation of A)")
    exp4_peaks = [50, 460]
    print(f"Observed peaks at {exp4_peaks[0]} kDa and {exp4_peaks[1]} kDa.")
    print(f"The {exp4_peaks[0]} kDa peak is the Protein A homodimer ({native_masses['A']} kDa). Dephosphorylation restores its native dimeric state.")
    print(f"The {exp4_peaks[1]} kDa peak shows the B-C-D complex ({native_masses['B'] + native_masses['C'] + native_masses['D']} kDa) remains intact.")
    print("Conclusion for Exp 4: The B-C-D complex is stable and non-phosphorylated Protein A cannot displace Protein B once it has bound.")
    print("-" * 60)

    # --- Step 5: Final Conclusion based on all data ---
    print("Step 5: Evaluating Answer Choices")
    print("\nBased on the analysis, let's evaluate statement G:")
    print("'Protein B never has a higher affinity for protein C or D than nonphosphorylated protein A, protein B always exists as a dimer in opposition to proteins A, C, and D.'")
    print("\n- Part 1: 'Protein B never has a higher affinity for protein C or D than nonphosphorylated protein A'")
    print("  This is SUPPORTED by Experiment 2, where the A(dimer)-C-D complex (210 kDa) formed in preference to a B-containing complex.")
    
    print("\n- Part 2: 'protein B always exists as a dimer in opposition to proteins A, C, and D.'")
    print("  This is SUPPORTED by the full dataset:")
    print(f"  - Protein B is always a dimer ({native_masses['B']} kDa).")
    print(f"  - Protein A exists as a dimer ({native_masses['A']} kDa) and a monomer ({protein_masses['A']} kDa).")
    print(f"  - Protein C is a monomer ({native_masses['C']} kDa).")
    print(f"  - Protein D is a monomer ({native_masses['D']} kDa).")
    print("  Therefore, Protein B's oligomeric state (always a dimer) is distinct from the others.")

    print("\nOverall Conclusion: Both parts of statement G are factually correct based on the experiments.")

if __name__ == '__main__':
    analyze_protein_interactions()