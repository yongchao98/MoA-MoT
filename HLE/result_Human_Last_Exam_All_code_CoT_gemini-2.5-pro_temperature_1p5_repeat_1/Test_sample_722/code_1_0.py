def solve_protein_puzzle():
    """
    Analyzes protein interaction data from SEC-MALS experiments to determine the correct conclusion.
    """
    # Theoretical masses (in kDa)
    A_theo = 25
    B_theo = 150
    C_theo = 60
    D_theo = 100
    Kinase = 40

    # --- Step 1: Analysis of Experiment 1 (Individual Proteins) ---
    print("--- Step 1: Analysis of Individual Proteins (Experiment 1) ---")
    A_exp1 = 50
    B_exp1 = 300
    C_exp1 = 60
    D_exp1 = 100

    # Determine oligomeric state of Protein A
    a_ratio = A_exp1 / A_theo
    print(f"Protein A: Theoretical mass is {A_theo} kDa. Experimental mass is {A_exp1} kDa.")
    print(f"The experimental mass ({A_exp1}) is {a_ratio:.0f} times the theoretical mass ({A_theo}). Conclusion: Protein A is a homodimer.")
    A_oligomer_mass = A_exp1

    # Determine oligomeric state of Protein B
    b_ratio = B_exp1 / B_theo
    print(f"Protein B: Theoretical mass is {B_theo} kDa. Experimental mass is {B_exp1} kDa.")
    print(f"The experimental mass ({B_exp1}) is {b_ratio:.0f} times the theoretical mass ({B_theo}). Conclusion: Protein B is a homodimer.")
    B_oligomer_mass = B_exp1

    # Determine oligomeric state of Protein C
    c_ratio = C_exp1 / C_theo
    print(f"Protein C: Theoretical mass is {C_theo} kDa. Experimental mass is {C_exp1} kDa.")
    print(f"The experimental mass ({C_exp1}) is {c_ratio:.0f} times the theoretical mass ({C_theo}). Conclusion: Protein C is a monomer.")
    C_oligomer_mass = C_exp1
    
    # Determine oligomeric state of Protein D
    d_ratio = D_exp1 / D_theo
    print(f"Protein D: Theoretical mass is {D_theo} kDa. Experimental mass is {D_exp1} kDa.")
    print(f"The experimental mass ({D_exp1}) is {d_ratio:.0f} times the theoretical mass ({D_theo}). Conclusion: Protein D is a monomer.")
    D_oligomer_mass = D_exp1
    print("-" * 50)

    # --- Step 2: Analysis of Experiment 2 (Protein Mixture) ---
    print("\n--- Step 2: Analysis of Protein Mixture (Experiment 2) ---")
    print("Observed peaks: 300 kDa and 210 kDa.")
    peak1_exp2 = 300
    peak2_exp2 = 210
    
    print(f"The peak at {peak1_exp2} kDa matches the mass of the Protein B dimer ({B_oligomer_mass} kDa), indicating it is unbound.")
    
    # Check the composition of the second peak
    complex_acd_mass = A_oligomer_mass + C_oligomer_mass + D_oligomer_mass
    print(f"The peak at {peak2_exp2} kDa matches the combined mass of the Protein A dimer, Protein C monomer, and Protein D monomer.")
    print(f"Final Equation: {A_oligomer_mass} (A dimer) + {C_oligomer_mass} (C monomer) + {D_oligomer_mass} (D monomer) = {complex_acd_mass} kDa")
    print("Conclusion: Unphosphorylated Protein A has a higher affinity for proteins C and D than Protein B does.")
    print("-" * 50)
    
    # --- Step 3: Analysis of Experiment 3 (Mixture + Kinase) ---
    print("\n--- Step 3: Analysis of Mixture + Kinase (Experiment 3) ---")
    print("Observed peaks: 25 kDa, 40 kDa, and 460 kDa.")
    peak1_exp3 = 25
    peak2_exp3 = 40
    peak3_exp3 = 460
    
    print(f"The peak at {peak1_exp3} kDa matches the mass of monomeric Protein A ({A_theo} kDa).")
    print(f"The peak at {peak2_exp3} kDa matches the mass of the Kinase ({Kinase} kDa), which is unbound.")
    
    # Check the composition of the third peak
    complex_bcd_mass = B_oligomer_mass + C_oligomer_mass + D_oligomer_mass
    print(f"The peak at {peak3_exp3} kDa matches the combined mass of the Protein B dimer, Protein C monomer, and Protein D monomer.")
    print(f"Final Equation: {B_oligomer_mass} (B dimer) + {C_oligomer_mass} (C monomer) + {D_oligomer_mass} (D monomer) = {complex_bcd_mass} kDa")
    print("Conclusion: Phosphorylation of Protein A causes it to dissociate into monomers and lose affinity for C and D, allowing Protein B to bind them instead.")
    print("-" * 50)

    # --- Step 4: Analysis of Experiment 4 (Dephosphorylation) ---
    print("\n--- Step 4: Analysis of Dephosphorylated Mixture (Experiment 4) ---")
    print("Observed peaks: 50 kDa and 460 kDa.")
    print(f"The peak at 50 kDa shows that dephosphorylated Protein A has re-formed its dimer ({A_oligomer_mass} kDa).")
    print(f"The peak at 460 kDa shows that the complex of B, C, and D ({complex_bcd_mass} kDa) remains stable.")
    print("Conclusion: Confirms Protein A was the target of phosphorylation and that it exists as a dimer in its unphosphorylated state.")
    print("-" * 50)

    # --- Step 5: Evaluate Answer Choices ---
    print("\n--- Step 5: Evaluating the Answer Choices ---")
    print("A: Incorrect. Phosphorylation targets Protein A, not B.")
    print("B: Incorrect. Phosphorylation of A *decreases* its affinity for C and D.")
    print("C: Incorrect. Phosphorylation targets Protein A, not B.")
    print("D: Incorrect. Phosphorylation of A *decreases* its affinity.")
    print("E: Incorrect. Protein A is a monomer when phosphorylated.")
    print("F: Incorrect. Protein A is a dimer when unphosphorylated. This statement is incomplete/misleading.")
    print("G: Incorrect. Protein A is also a dimer, so B is not in 'opposition' to A's oligomeric state.")
    print("H: Incorrect. Phosphorylation targets Protein A, not B.")
    print("I: Incorrect. C and H are both false.")
    print("\nFinal Verdict: Every answer from A to I contains a statement that contradicts the experimental data.")

if __name__ == '__main__':
    solve_protein_puzzle()
<<<J>>>