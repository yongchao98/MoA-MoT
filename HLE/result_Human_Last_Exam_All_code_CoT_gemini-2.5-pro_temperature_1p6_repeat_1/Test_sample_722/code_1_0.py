def analyze_protein_interactions():
    """
    Analyzes SEC-MALS data to determine protein interactions and select the correct conclusion.
    """
    # Define initial protein masses
    protein_mass_theoretical = {
        'A': 25,
        'B': 150,
        'C': 60,
        'D': 100,
    }
    kinase_mass = 40

    print("Analysis of Protein Interaction Experiments")
    print("==========================================")
    print()

    # --- Step 1: Analyze Experiment 1 ---
    print("--- Analysis of Experiment 1: Individual Proteins ---")
    exp1_masses = {'A': 50, 'B': 300, 'C': 60, 'D': 100}
    A_oligomer = exp1_masses['A'] / protein_mass_theoretical['A']
    B_oligomer = exp1_masses['B'] / protein_mass_theoretical['B']
    C_oligomer = exp1_masses['C'] / protein_mass_theoretical['C']
    D_oligomer = exp1_masses['D'] / protein_mass_theoretical['D']

    print(f"Protein A: Theoretical Mass = {protein_mass_theoretical['A']} kDa, Observed Mass = {exp1_masses['A']} kDa.")
    print(f"Calculation: {exp1_masses['A']} / {protein_mass_theoretical['A']} = {int(A_oligomer)}. Conclusion: Protein A is a homodimer.")
    A_dimer_mass = exp1_masses['A']
    print()

    print(f"Protein B: Theoretical Mass = {protein_mass_theoretical['B']} kDa, Observed Mass = {exp1_masses['B']} kDa.")
    print(f"Calculation: {exp1_masses['B']} / {protein_mass_theoretical['B']} = {int(B_oligomer)}. Conclusion: Protein B is a homodimer.")
    B_dimer_mass = exp1_masses['B']
    print()

    print(f"Protein C: Theoretical Mass = {protein_mass_theoretical['C']} kDa, Observed Mass = {exp1_masses['C']} kDa.")
    print(f"Calculation: {exp1_masses['C']} / {protein_mass_theoretical['C']} = {int(C_oligomer)}. Conclusion: Protein C is a monomer.")
    C_mono_mass = exp1_masses['C']
    print()

    print(f"Protein D: Theoretical Mass = {protein_mass_theoretical['D']} kDa, Observed Mass = {exp1_masses['D']} kDa.")
    print(f"Calculation: {exp1_masses['D']} / {protein_mass_theoretical['D']} = {int(D_oligomer)}. Conclusion: Protein D is a monomer.")
    D_mono_mass = exp1_masses['D']
    print("\n")


    # --- Step 2: Analyze Experiment 2 ---
    print("--- Analysis of Experiment 2: Mixture of A, B, C, D ---")
    exp2_peak1 = 300
    exp2_peak2 = 210
    print(f"Observed peaks at {exp2_peak1} kDa and {exp2_peak2} kDa.")
    print(f"The {exp2_peak1} kDa peak corresponds to the mass of the Protein B dimer ({B_dimer_mass} kDa), which is free in solution.")
    print("Let's test if the second peak is a complex of the remaining proteins: A (dimer), C (monomer), and D (monomer).")
    complex1_mass_calc = A_dimer_mass + C_mono_mass + D_mono_mass
    print(f"Calculation: Mass(A_dimer) + Mass(C_mono) + Mass(D_mono) = {A_dimer_mass} + {C_mono_mass} + {D_mono_mass} = {complex1_mass_calc} kDa.")
    print(f"This calculated mass ({complex1_mass_calc} kDa) matches the second peak of {exp2_peak2} kDa.")
    print("Conclusion: When mixed, a complex of A-C-D forms, leaving B free. This indicates that dephosphorylated Protein A (dimer) has a higher affinity for proteins C and D than Protein B (dimer) does.")
    print("\n")

    # --- Step 3: Analyze Experiment 3 ---
    print("--- Analysis of Experiment 3: Mixture + Kinase ---")
    exp3_peak1 = 25
    exp3_peak2 = 40
    exp3_peak3 = 460
    print(f"Observed peaks at {exp3_peak1} kDa, {exp3_peak2} kDa, and {exp3_peak3} kDa.")
    print(f"The {exp3_peak2} kDa peak matches the mass of the kinase ({kinase_mass} kDa). The kinase is free.")
    print(f"The {exp3_peak1} kDa peak matches the theoretical mass of Protein A monomer ({protein_mass_theoretical['A']} kDa).")
    print("This implies the kinase phosphorylated Protein A, causing its dimer to dissociate into monomers.")
    print(f"Let's test if the large {exp3_peak3} kDa peak is a complex of B (dimer), C (monomer), and D (monomer).")
    complex2_mass_calc = B_dimer_mass + C_mono_mass + D_mono_mass
    print(f"Calculation: Mass(B_dimer) + Mass(C_mono) + Mass(D_mono) = {B_dimer_mass} + {C_mono_mass} + {D_mono_mass} = {complex2_mass_calc} kDa.")
    print(f"This calculated mass ({complex2_mass_calc} kDa) matches the third peak of {exp3_peak3} kDa.")
    print("Conclusion: Phosphorylation of Protein A causes it to become a monomer and prevents it from binding to C and D. This allows Protein B to form a complex with C and D.")
    print("\n")

    # --- Step 4: Analyze Experiment 4 ---
    print("--- Analysis of Experiment 4: After dephosphorylation of A ---")
    exp4_peak1 = 50
    exp4_peak2 = 460
    print(f"Observed peaks at {exp4_peak1} kDa and {exp4_peak2} kDa.")
    print(f"The {exp4_peak1} kDa peak is the Protein A dimer ({A_dimer_mass} kDa), which has re-formed after dephosphorylation.")
    print(f"The {exp4_peak2} kDa peak is the stable B-C-D complex ({complex2_mass_calc} kDa) from Experiment 3.")
    print("Conclusion: The B-C-D complex, once formed, is kinetically stable and is not displaced by dephosphorylated Protein A, despite A's higher intrinsic affinity for C+D.")
    print("\n")

    # --- Step 5: Evaluate Answer Choices ---
    print("--- Final Evaluation of Answer Choices ---")
    print("Based on the experimental analysis:")
    print("A: Incorrect. No evidence protein B is phosphorylated.")
    print("B: Incorrect. Phosphorylation of A *decreases* its affinity for the C+D complex.")
    print("C: Incorrect. No evidence protein B is phosphorylated.")
    print("D: Incorrect. Phosphorylation of A *decreases* affinity.")
    print("E: Incorrect. Protein A is a monomer when phosphorylated.")
    print("F: Incorrect. Protein A also exists as a dimer.")
    print("G: Correct. Experiment 2 shows nonphosphorylated protein A has a higher affinity for C+D than B. Protein B is always a dimer, unlike A (monomer/dimer) and C, D (monomers).")
    print("H: Incorrect. No evidence protein B is phosphorylated.")
    print("I: Incorrect, as C and H are incorrect.")
    print("J: Incorrect, as G is correct.")
    print()

    final_answer = 'G'
    print(f'<<<{final_answer}>>>')

# Run the analysis
analyze_protein_interactions()