def analyze_protein_interactions():
    """
    Analyzes the provided experimental data to determine protein interactions and evaluates the given answer choices.
    """

    # --- Step 1: Define Masses and Analyze Baseline from Experiment 1 ---
    print("--- Analysis Step 1: Determining Protein States from Experiment 1 ---")
    protein_A_mono = 25
    protein_B_mono = 150
    protein_C_mono = 60
    protein_D_mono = 100
    kinase_mass = 40

    # From Exp 1 Results
    protein_A_exp1 = 50
    protein_B_exp1 = 300
    protein_C_exp1 = 60
    protein_D_exp1 = 100

    A_dimer = protein_A_mono * 2
    B_dimer = protein_B_mono * 2

    print(f"Protein A: Theoretical monomer is {protein_A_mono} kDa. Observed is {protein_A_exp1} kDa.")
    print(f"Calculation: {protein_A_mono} * 2 = {A_dimer}. Conclusion: Protein A is a homodimer.")
    print(f"Protein B: Theoretical monomer is {protein_B_mono} kDa. Observed is {protein_B_exp1} kDa.")
    print(f"Calculation: {protein_B_mono} * 2 = {B_dimer}. Conclusion: Protein B is a homodimer.")
    print(f"Protein C: Theoretical monomer {protein_C_mono} kDa, observed {protein_C_exp1} kDa. Conclusion: Protein C is a monomer.")
    print(f"Protein D: Theoretical monomer {protein_D_mono} kDa, observed {protein_D_exp1} kDa. Conclusion: Protein D is a monomer.\n")

    # --- Step 2: Analyze Experiment 2 ---
    print("--- Analysis Step 2: Interpreting Protein Complex in Experiment 2 ---")
    peak1_exp2 = 300
    peak2_exp2 = 210

    print(f"Observed peaks are {peak1_exp2} kDa and {peak2_exp2} kDa.")
    print(f"The {peak1_exp2} kDa peak corresponds to the free Protein B dimer ({B_dimer} kDa).")
    
    # Check hypothesis for the 210 kDa peak
    complex_A_C_D = A_dimer + protein_C_mono + protein_D_mono
    print(f"Hypothesis for {peak2_exp2} kDa peak: Protein A dimer + Protein C + Protein D.")
    print(f"Calculation: {A_dimer} + {protein_C_mono} + {protein_D_mono} = {complex_A_C_D} kDa.")
    print(f"Conclusion: The 210 kDa peak is a complex of A-dimer, C, and D. This means non-phosphorylated A has a higher affinity for C+D than B does.\n")

    # --- Step 3: Analyze Experiment 3 ---
    print("--- Analysis Step 3: Interpreting the Effect of a Kinase in Experiment 3 ---")
    peak1_exp3 = 25
    peak2_exp3 = 40
    peak3_exp3 = 460
    print(f"Observed peaks are {peak1_exp3} kDa, {peak2_exp3} kDa, and {peak3_exp3} kDa.")
    print(f"The {peak1_exp3} kDa peak corresponds to the Protein A monomer ({protein_A_mono} kDa). This indicates phosphorylation caused the A-dimer to dissociate.")
    print(f"The {peak2_exp3} kDa peak corresponds to the free kinase ({kinase_mass} kDa).")

    # Check hypothesis for the 460 kDa peak
    complex_B_C_D = B_dimer + protein_C_mono + protein_D_mono
    print(f"Hypothesis for {peak3_exp3} kDa peak: Protein B dimer + Protein C + Protein D.")
    print(f"Calculation: {B_dimer} + {protein_C_mono} + {protein_D_mono} = {complex_B_C_D} kDa.")
    print("Conclusion: Phosphorylation of A decreases its affinity for C+D, allowing B to form a complex with them.\n")

    # --- Step 4: Analyze Experiment 4 ---
    print("--- Analysis Step 4: Interpreting Dephosphorylation in Experiment 4 ---")
    peak1_exp4 = 50
    peak2_exp4 = 460
    print(f"Observed peaks are {peak1_exp4} kDa and {peak2_exp4} kDa.")
    print(f"The sample started with the stable B-C-D complex from Exp 3 and dephosphorylated A-monomers.")
    print(f"The {peak1_exp4} kDa peak is the re-formed Protein A dimer ({A_dimer} kDa).")
    print(f"The {peak2_exp4} kDa peak shows the B-C-D complex ({complex_B_C_D} kDa) remained stable.")
    print("Conclusion: The B-C-D complex, once formed, is very stable and cannot be displaced by non-phosphorylated Protein A. Therefore, non-phosphorylated A does NOT always have a higher affinity.\n")

    # --- Step 5: Evaluate Answer Choices ---
    print("--- Step 5: Evaluating Answer Choices ---")
    print("A: Incorrect. Protein A is phosphorylated, not B.")
    print("B: Incorrect. Phosphorylation of A decreases its affinity for C and D.")
    print("C: Incorrect. Protein A is phosphorylated, not B, and non-phosphorylated A does not always have higher affinity (see Exp 4).")
    print("D: Incorrect. Phosphorylation of A decreases its affinity.")
    print("E: Incorrect. Protein A becomes a monomer when phosphorylated (Exp 3).")
    print("F: Incorrect. Protein A is usually a dimer, not a monomer.")
    print("G: Incorrect. The persistence of the B-C-D complex in Exp 4 shows non-phosphorylated A does not always win the competition.")
    print("H: Incorrect. Protein A is phosphorylated, not B.")
    print("I: Incorrect, as C and H are incorrect.")
    print("J: Correct. All other statements contain flaws based on the experimental evidence.")

if __name__ == "__main__":
    analyze_protein_interactions()
    print("\n<<<J>>>")