import sys

def solve_protein_puzzle():
    """
    Analyzes protein interaction data from SEC-MALS experiments to determine the correct conclusion.
    """
    # --- Step 1: Define theoretical masses and analyze Experiment 1 ---
    print("--- Analysis of Experiment 1: Determining Native Protein States ---")
    protein_A_mono = 25
    protein_B_mono = 150
    protein_C_mono = 60
    protein_D_mono = 100
    kinase_mass = 40

    # Observed masses from Experiment 1
    obs_A = 50
    obs_B = 300
    obs_C = 60
    obs_D = 100

    print(f"Theoretical masses (kDa): A={protein_A_mono}, B={protein_B_mono}, C={protein_C_mono}, D={protein_D_mono}")
    print(f"Observed masses (kDa): A={obs_A}, B={obs_B}, C={obs_C}, D={obs_D}\n")

    # Determine oligomeric states
    state_A = obs_A / protein_A_mono
    state_B = obs_B / protein_B_mono
    
    print(f"Result for Protein A: Observed mass is {obs_A} kDa. The calculation is {obs_A} / {protein_A_mono} = {state_A}. Protein A is a homodimer (A2).")
    print(f"Result for Protein B: Observed mass is {obs_B} kDa. The calculation is {obs_B} / {protein_B_mono} = {state_B}. Protein B is a homodimer (B2).")
    print(f"Result for Protein C: Observed mass ({obs_C} kDa) equals theoretical mass. Protein C is a monomer.")
    print(f"Result for Protein D: Observed mass ({obs_D} kDa) equals theoretical mass. Protein D is a monomer.")
    print("-" * 50)

    # --- Step 2: Analyze Experiment 2 ---
    print("\n--- Analysis of Experiment 2: Mixture of Native Proteins ---")
    peak1_exp2 = 300
    peak2_exp2 = 210
    print(f"Observed peaks at {peak1_exp2} kDa and {peak2_exp2} kDa.")
    
    # Identify peaks
    complex_A2_C_D = obs_A + obs_C + obs_D
    print(f"The {peak1_exp2} kDa peak matches the mass of the Protein B homodimer (B2), which is unbound.")
    print(f"The {peak2_exp2} kDa peak matches the mass of a complex of Protein A dimer + Protein C + Protein D. The calculation is {obs_A} + {obs_C} + {obs_D} = {complex_A2_C_D} kDa.")
    print("Conclusion: Non-phosphorylated Protein A (as A2 dimer) has a higher affinity for proteins C and D than the Protein B dimer (B2).")
    print("-" * 50)

    # --- Step 3: Analyze Experiment 3 ---
    print("\n--- Analysis of Experiment 3: Mixture with Kinase (Phosphorylation) ---")
    peak1_exp3 = 25
    peak2_exp3 = 40
    peak3_exp3 = 460
    print(f"Observed peaks at {peak1_exp3} kDa, {peak2_exp3} kDa, and {peak3_exp3} kDa.")
    
    # Identify peaks
    complex_B2_C_D = obs_B + obs_C + obs_D
    print(f"The {peak1_exp3} kDa peak matches the theoretical mass of monomeric Protein A ({protein_A_mono} kDa). This indicates phosphorylation breaks the A2 dimer.")
    print(f"The {peak2_exp3} kDa peak matches the mass of the free kinase ({kinase_mass} kDa).")
    print(f"The {peak3_exp3} kDa peak matches a complex of Protein B dimer + Protein C + Protein D. The calculation is {obs_B} + {obs_C} + {obs_D} = {complex_B2_C_D} kDa.")
    print("Conclusion: When Protein A is phosphorylated, it becomes a monomer and loses its affinity for C and D, allowing the B2 dimer to bind them instead.")
    print("-" * 50)

    # --- Step 4 & 5: Analyze Experiment 4 and Evaluate Answer Choice G ---
    print("\n--- Final Evaluation based on all Experiments ---")
    print("Answer choice G states: 'Protein B never has a higher affinity for protein C or D than nonphosphorylated protein A, protein B always exists as a dimer in opposition to proteins A, C, and D.'")
    
    # Evaluate first part of G
    print("\nPart 1 of G: 'Protein B never has a higher affinity for protein C or D than nonphosphorylated protein A'")
    print("This is proven by Experiment 2. When all native (nonphosphorylated) proteins were mixed, the A2-C-D complex formed, leaving B2 free. This shows A2's affinity for C+D is higher.")
    
    # Evaluate second part of G
    print("\nPart 2 of G: 'protein B always exists as a dimer in opposition to proteins A, C, and D.'")
    print(f"This is also correct. B is consistently a dimer ({obs_B} kDa). In contrast, A can be a dimer ({obs_A} kDa) or a monomer ({protein_A_mono} kDa), while C and D are always monomers.")
    
    print("\nOverall Conclusion: Both parts of statement G are correct and supported by the experimental evidence.")

if __name__ == "__main__":
    solve_protein_puzzle()
    # The final answer is derived from the logical conclusions presented above.
    sys.stdout.flush() # Ensure all prints are displayed before the final answer tag.
    print("\n<<<G>>>")