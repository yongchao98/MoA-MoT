def analyze_protein_interactions():
    """
    Analyzes protein interaction data from SEC-MALS experiments and prints the reasoning.
    """
    # Theoretical masses
    mass_A_mono = 25
    mass_B_mono = 150
    mass_C_mono = 60
    mass_D_mono = 100
    mass_kinase = 40

    print("--- Analysis of Protein-Protein Interactions ---")

    print("\n[Experiment 1: Individual Proteins]")
    mass_A_exp1 = 50
    mass_B_exp1 = 300
    print(f"Protein A: Measured mass is {mass_A_exp1} kDa. This corresponds to a dimer.")
    print(f"Equation: {mass_A_mono} kDa * 2 = {mass_A_exp1} kDa")
    print(f"Protein B: Measured mass is {mass_B_exp1} kDa. This corresponds to a dimer.")
    print(f"Equation: {mass_B_mono} kDa * 2 = {mass_B_exp1} kDa")
    print(f"Protein C: Measured mass is {mass_C_mono} kDa, a monomer.")
    print(f"Protein D: Measured mass is {mass_D_mono} kDa, a monomer.")
    
    # Inferred oligomeric masses for complexes
    mass_A_dimer = mass_A_exp1
    mass_B_dimer = mass_B_exp1

    print("\n[Experiment 2: All Proteins Mixed (Non-phosphorylated)]")
    mass_complex_A_C_D = 210
    print(f"Peak 1: 300 kDa (Free Protein B dimer).")
    print(f"Peak 2: {mass_complex_A_C_D} kDa. This is the A-C-D complex.")
    print(f"Equation for Peak 2: {mass_A_dimer} (A_dimer) + {mass_C_mono} (C) + {mass_D_mono} (D) = {mass_A_dimer + mass_C_mono + mass_D_mono} kDa")
    print("Conclusion: Dephosphorylated Protein A (dimer) has a higher affinity for C+D than Protein B (dimer).")

    print("\n[Experiment 3: All Proteins + Kinase (Phosphorylation of A)]")
    mass_complex_B_C_D = 460
    print(f"Peak 1: {mass_A_mono} kDa (Phosphorylated Protein A monomer).")
    print(f"Peak 2: {mass_kinase} kDa (Free Kinase).")
    print(f"Peak 3: {mass_complex_B_C_D} kDa. This is the B-C-D complex.")
    print(f"Equation for Peak 3: {mass_B_dimer} (B_dimer) + {mass_C_mono} (C) + {mass_D_mono} (D) = {mass_B_dimer + mass_C_mono + mass_D_mono} kDa")
    print("Conclusion: Phosphorylation causes Protein A to dissociate into monomers and lose its affinity for C+D.")

    print("\n[Experiment 4: Dephosphorylation of A]")
    print(f"Peak 1: {mass_A_dimer} kDa (Dephosphorylated Protein A re-forms its dimer).")
    print(f"Peak 2: {mass_complex_B_C_D} kDa (The B-C-D complex remains stable).")
    print(f"Equation for Peak 2: {mass_B_dimer} (B_dimer) + {mass_C_mono} (C) + {mass_D_mono} (D) = {mass_complex_B_C_D} kDa")
    print("Conclusion: Dephosphorylated A cannot displace B from the pre-formed complex.")
    
    print("\n--- Final Conclusion ---")
    print("The evidence supports option G:")
    print("1. Protein B never has higher affinity for C/D than nonphosphorylated A (Exp 2).")
    print("2. Protein B always exists as a dimer (All Exps).")
    print("3. B's dimer state is in opposition to A (monomer/dimer), C (monomer), and D (monomer).")


# Execute the analysis
analyze_protein_interactions()

print("<<<G>>>")