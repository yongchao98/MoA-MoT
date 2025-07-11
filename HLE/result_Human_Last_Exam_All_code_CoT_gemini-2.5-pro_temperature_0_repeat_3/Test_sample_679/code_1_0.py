def solve_nmr_puzzle():
    """
    Analyzes the provided molecular formula and 13C NMR data to determine the IUPAC name of the compound.
    """
    
    # Given data
    molecular_formula = "C7H14"
    nmr_signals = {
        145: 's', 112: 't', 48: 't', 27: 'd', 22: 'q', 21: 'q'
    }
    
    print("Step-by-step analysis to determine the IUPAC name:")
    print("-" * 50)
    
    # Step 1: Analyze Molecular Formula
    print("1. Molecular Formula Analysis:")
    print(f"   - The formula is {molecular_formula}. The general formula for a saturated alkane is CnH(2n+2), which would be C7H16.")
    print("   - The compound has two fewer hydrogens, indicating one degree of unsaturation (one double bond or one ring).\n")
    
    # Step 2: Analyze NMR Data
    print("2. 13C NMR Data Analysis:")
    print("   - The signals at 145(s) and 112(t) are in the alkene region (>100 ppm).")
    print("   - 145(s) corresponds to a quaternary carbon (C) with no hydrogens.")
    print("   - 112(t) corresponds to a methylene carbon (CH2) with two hydrogens.")
    print("   - This confirms the presence of a terminal alkene group: >C=CH2.\n")
    
    # Step 3: Assemble Fragments
    print("3. Assembling the Structure:")
    print("   - The molecule has 7 carbons but only 6 NMR signals, which means two carbons are chemically equivalent.")
    print("   - The remaining signals are 48(t), 27(d), 22(q), and 21(q).")
    print("   - These correspond to an ethyl group (-CH2CH3) and an isopropyl group (-CH(CH3)2).")
    print("   - The two methyl groups of the isopropyl group are equivalent, accounting for the 7th carbon.\n")
    
    # Step 4: Final Structure and Name
    print("4. Final Structure and IUPAC Name:")
    print("   - Attaching the ethyl and isopropyl groups to the quaternary carbon of the alkene gives the following structure:")
    print("""
          CH3
           |
      CH3-CH--C==CH2
              |
              CH2-CH3
    """)
    print("   - The longest chain containing the double bond is a butene chain.")
    print("   - Numbering from the double bond, there is an ethyl group at C2 and a methyl group at C3.")
    print("   - The IUPAC name is 2-ethyl-3-methylbut-1-ene.\n")

    # Step 5: Signal Assignment
    print("5. Signal Assignment:")
    print("   - 145(s): Quaternary alkene carbon (=C<)")
    print("   - 112(t): Terminal alkene carbon (=CH2)")
    print("   - 48(t): Methylene carbon of the ethyl group (-CH2-)")
    print("   - 27(d): Methine carbon of the isopropyl group (-CH-)")
    print("   - 22(q) & 21(q): The two signals for the methyl groups (-CH3)")
    print("-" * 50)
    
    final_answer = "2-ethyl-3-methylbut-1-ene"
    print(f"The IUPAC name of the compound is: {final_answer}")

solve_nmr_puzzle()