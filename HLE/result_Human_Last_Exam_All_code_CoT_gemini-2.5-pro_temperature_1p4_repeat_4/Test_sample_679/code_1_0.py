def solve_structure():
    """
    This function analyzes the provided molecular formula and 13C NMR data
    to determine the IUPAC name of the hydrocarbon.
    """

    # Step 1: Analyze molecular formula and Degree of Unsaturation (DBE).
    # Formula: C7H14 -> DBE = 1 (one double bond or one ring).
    molecular_formula = "C7H14"
    dbe = 1

    # Step 2: Analyze 13C NMR data.
    # Signals: 145(s), 112(t), 48(t), 27(d), 22(q), 21(q).
    # 6 signals for 7 carbons means two carbons are chemically equivalent.
    # Signals at 145 and 112 ppm confirm a C=C double bond.
    nmr_signals = {
        145: 's', 112: 't', 48: 't',
        27: 'd', 22: 'q', 21: 'q'
    }

    # Step 3: Piece together the fragments based on multiplicities and chemical shifts.
    # 145(s) and 112(t) in the alkene region -> >C=CH2 fragment.
    # The remaining 5 carbons and 12 hydrogens form two alkyl groups attached to the quaternary carbon.
    # 27(d) and 21(q, for 2C) -> isopropyl group [-CH(CH3)2].
    # 48(t) and 22(q) -> ethyl group [-CH2CH3].
    fragment1 = "ethyl group"
    fragment2 = "isopropyl group"
    backbone = ">C=CH2"

    # Step 4: Assemble the final structure and determine the IUPAC name.
    # Structure: The ethyl and isopropyl groups are attached to the quaternary carbon of the >C=CH2 fragment.
    # IUPAC Naming:
    #   - Longest chain with C=C: butene.
    #   - Double bond at position 1: but-1-ene.
    #   - Substituents: ethyl at C2, methyl at C3.
    # Final Name: 2-ethyl-3-methylbut-1-ene.

    final_iupac_name = "2-ethyl-3-methylbut-1-ene"
    
    # The problem asks to output the numbers in the final name/equation.
    # The numbers in the name "2-ethyl-3-methylbut-1-ene" are 2, 3, and 1.
    # Printing the full name will display these numbers.
    print(f"The analysis of the molecular formula and 13C NMR data leads to the following structure.")
    print(f"Structure fragments: A >C=CH2 group, an ethyl group, and an isopropyl group.")
    print(f"Final Assembled Structure: CH2=C(CH2CH3)(CH(CH3)2)")
    print(f"\nThe IUPAC name of the compound is:")
    print(final_iupac_name)

solve_structure()
<<<2-ethyl-3-methylbut-1-ene>>>