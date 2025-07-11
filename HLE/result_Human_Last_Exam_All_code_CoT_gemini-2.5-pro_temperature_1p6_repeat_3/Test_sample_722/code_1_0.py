def analyze_protein_interactions():
    """
    Analyzes the SEC-MALS data step-by-step and prints the logical deductions.
    """
    # Theoretical monomer masses in kDa
    mass_A_mono = 25
    mass_B_mono = 150
    mass_C_mono = 60
    mass_D_mono = 100
    mass_kinase = 40

    print("--- Step 1: Analysis of Experiment 1 (Individual Proteins) ---")
    # Observed masses in Experiment 1
    obs_A = 50
    obs_B = 300
    obs_C = 60
    obs_D = 100

    # Determine oligomeric states
    state_A = obs_A / mass_A_mono
    state_B = obs_B / mass_B_mono
    state_C = obs_C / mass_C_mono
    state_D = obs_D / mass_D_mono

    print(f"Protein A: Observed {obs_A} kDa / Theoretical {mass_A_mono} kDa = {int(state_A)}. Conclusion: Protein A is a Homodimer (A₂).")
    print(f"Protein B: Observed {obs_B} kDa / Theoretical {mass_B_mono} kDa = {int(state_B)}. Conclusion: Protein B is a Homodimer (B₂).")
    print(f"Protein C: Observed {obs_C} kDa / Theoretical {mass_C_mono} kDa = {int(state_C)}. Conclusion: Protein C is a Monomer.")
    print(f"Protein D: Observed {obs_D} kDa / Theoretical {mass_D_mono} kDa = {int(state_D)}. Conclusion: Protein D is a Monomer.")
    print("-" * 60)

    # Use the determined oligomeric masses for subsequent calculations
    mass_A_oligomer = obs_A
    mass_B_oligomer = obs_B

    print("\n--- Step 2: Analysis of Experiment 2 (Mixture of A, B, C, D) ---")
    # Observed peaks in Experiment 2
    peak1_exp2 = 300
    peak2_exp2 = 210
    print(f"Observed Peaks: {peak1_exp2} kDa and {peak2_exp2} kDa.")
    print(f"The {peak1_exp2} kDa peak corresponds to the free Protein B dimer ({mass_B_oligomer} kDa).")
    print("Let's check if the remaining proteins form the second peak:")
    complex_ACD_mass = mass_A_oligomer + mass_C_mono + mass_D_mono
    print(f"Equation: Mass(A dimer) + Mass(C monomer) + Mass(D monomer) = {mass_A_oligomer} + {mass_C_mono} + {mass_D_mono} = {complex_ACD_mass} kDa.")
    print(f"This calculated mass ({complex_ACD_mass} kDa) matches the second peak ({peak2_exp2} kDa).")
    print("Conclusion: A₂CD complex forms, leaving B₂ free. This means non-phosphorylated Protein A (as a dimer) has a higher affinity for C+D than Protein B does.")
    print("-" * 60)

    print("\n--- Step 3: Analysis of Experiment 3 (Mixture + Kinase) ---")
    # Observed peaks in Experiment 3
    peak1_exp3 = 25
    peak2_exp3 = 40
    peak3_exp3 = 460
    print(f"Observed Peaks: {peak1_exp3} kDa, {peak2_exp3} kDa, and {peak3_exp3} kDa.")
    print(f"The {peak2_exp3} kDa peak corresponds to the free Kinase ({mass_kinase} kDa).")
    print(f"The {peak1_exp3} kDa peak corresponds to monomeric Protein A ({mass_A_mono} kDa). This indicates phosphorylation of A causes its dimer to dissociate.")
    print("Let's check if the remaining proteins form the large complex:")
    complex_BCD_mass = mass_B_oligomer + mass_C_mono + mass_D_mono
    print(f"Equation: Mass(B dimer) + Mass(C monomer) + Mass(D monomer) = {mass_B_oligomer} + {mass_C_mono} + {mass_D_mono} = {complex_BCD_mass} kDa.")
    print(f"This calculated mass ({complex_BCD_mass} kDa) matches the third peak ({peak3_exp3} kDa).")
    print("Conclusion: Phosphorylation of Protein A decreases its affinity for C+D, allowing Protein B (dimer) to bind and form a B₂CD complex.")
    print("-" * 60)

    print("\n--- Final Evaluation of Answer Choices ---")
    print("Based on the analysis:")
    print("1. Non-phosphorylated Protein A (as A₂) has a higher affinity for C+D than Protein B (as B₂) does (from Exp 2).")
    print("2. Phosphorylation of Protein A causes it to become a monomer and decreases its affinity for C+D (from Exp 3).")
    print("3. Protein B is consistently a dimer (300 kDa), while Protein A's state changes (monomer/dimer) and C and D are monomers.")
    print("\nThis evidence strongly supports option G:")
    print("G. Protein B never has a higher affinity for protein C or D than nonphosphorylated protein A, protein B always exists as a dimer in opposition to proteins A, C, and D.")

analyze_protein_interactions()
<<<G>>>