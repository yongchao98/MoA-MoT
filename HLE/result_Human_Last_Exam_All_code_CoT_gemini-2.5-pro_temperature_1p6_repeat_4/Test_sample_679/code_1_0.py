def solve_structure():
    """
    This function solves for the IUPAC name of a hydrocarbon based on its
    molecular formula and 13C NMR data.
    """
    # Given Data
    molecular_formula = "C7H14"
    nmr_shifts = {
        145: 's', # singlet
        112: 't', # triplet
        48: 't',  # triplet
        27: 'd',  # doublet
        22: 'q',  # quartet
        21: 'q'   # quartet
    }

    print("Step 1: Analyze Molecular Formula")
    print(f"The formula is {molecular_formula}. This matches the general formula CnH2n, indicating one degree of unsaturation (one double bond or one ring).")
    print("-" * 20)

    print("Step 2: Analyze 13C NMR Spectrum")
    print(f"The given shifts are: {list(nmr_shifts.keys())}")
    print("The shifts at 145 ppm (s) and 112 ppm (t) are in the alkene region (100-150 ppm), confirming a C=C double bond.")
    print(" - The 145(s) signal indicates a quaternary carbon (C) in the double bond.")
    print(" - The 112(t) signal indicates a terminal CH2 carbon (=CH2) in the double bond.")
    print(" - This establishes the alkene structure as R2C=CH2.")
    print("There are 7 carbons in the formula but only 6 NMR signals. This means two carbons are chemically equivalent and share one signal.")
    print("-" * 20)
    
    print("Step 3: Deduce and Assemble Structure")
    print("Let's list the carbon types from the multiplicities:")
    print(" - 145(s): 1 x C (quaternary)")
    print(" - 112(t): 1 x CH2")
    print(" - 48(t):  1 x CH2")
    print(" - 27(d):  1 x CH")
    print(" - 22(q):  1 or 2 x CH3")
    print(" - 21(q):  1 or 2 x CH3")
    print("Counting these up gives 6 carbons. To get to 7, one signal must represent two carbons. Based on the commonality of equivalent methyl groups, we hypothesize one of the 'q' signals is for two equivalent CH3 groups.")
    print("The structure that fits all data points is 2-ethyl-3-methylbut-1-ene.")
    print("""
    Structure verification:
                CH2-CH3 (ethyl group)
               /
    H2C = C --- CH --- CH3
    1    2     3      4
               |
               CH3
    """)
    print("-" * 20)

    print("Step 4: Verify assignments and determine IUPAC Name")
    print("Let's assign the signals to the proposed structure:")
    print(f" - C1 (=CH2): A terminal alkene CH2. Signal: 112(t). Correct.")
    print(f" - C2 (=C<): A quaternary alkene carbon. Signal: 145(s). Correct.")
    print(f" - Ethyl -CH2-: Allylic to the double bond, so it's deshielded. Signal: 48(t). Correct.")
    print(f" - C3 (-CH<): A methine carbon. Signal: 27(d). Correct.")
    print(f" - Two CH3 groups at C3 and C4 are equivalent. Signal: 22(q) (represents 2 carbons). Correct.")
    print(f" - Ethyl -CH3: A methyl group. Signal: 21(q). Correct.")
    print("The proposed structure fits all the data.")
    print("-" * 20)

    final_name = "2-ethyl-3-methylbut-1-ene"
    print("Final Answer: The IUPAC name of the compound is:")
    print(final_name)

solve_structure()