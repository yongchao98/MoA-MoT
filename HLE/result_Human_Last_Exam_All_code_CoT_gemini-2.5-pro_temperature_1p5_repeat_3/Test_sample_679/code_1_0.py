def solve_nmr_puzzle():
    """
    This script explains the step-by-step reasoning to identify the compound
    with molecular formula C7H14 from its 13C NMR data.
    """
    molecular_formula = "C7H14"
    nmr_data = {
        145: "s",
        112: "t",
        48: "t",
        27: "d",
        22: "q",
        21: "q"
    }

    print("--- Step-by-Step Analysis ---")

    # Step 1: Analyze molecular formula
    print("\nStep 1: Analyze Molecular Formula")
    print(f"The molecular formula is {molecular_formula}.")
    print("This corresponds to a general formula of CnH2n, indicating one degree of unsaturation (a double bond or a ring).")

    # Step 2: Analyze 13C NMR Data
    print("\nStep 2: Analyze 13C NMR Data")
    print("The 13C NMR data shows 6 signals:")
    print("145 ppm (s): A quaternary sp2 carbon (>C=).")
    print("112 ppm (t): A methylene sp2 carbon (=CH2).")
    print("The signals at 145(s) and 112(t) confirm the presence of a C=C double bond, specifically a >C=CH2 group.")
    print("The remaining signals are from sp3 carbons:")
    print("48 ppm (t): A methylene carbon (-CH2-).")
    print("27 ppm (d): A methine carbon (-CH-).")
    print("22 ppm (q): A methyl carbon (-CH3).")
    print("21 ppm (q): Another methyl carbon (-CH3).")

    # Step 3: Propose and verify structure
    print("\nStep 3: Propose and Verify Structure")
    print("The molecule has 7 carbons but only 6 signals, which means two carbons are represented by a single peak.")
    print("The structure that fits all these fragments is 2-ethyl-3-methyl-1-butene.")
    print("This structure should theoretically have 7 unique carbon signals (1s, 2t, 1d, 3q).")
    print("The observed spectrum has 6 signals (1s, 2t, 1d, 2q).")

    # Step 4: Reconcile with data
    print("\nStep 4: Reconcile Prediction with Data")
    print("The discrepancy is resolved by the hint about signal overlapping.")
    print("The three methyl groups in the structure are chemically distinct, but two of them (the ones on the chiral center) are in very similar environments.")
    print("It is highly likely that two of the three methyl signals overlap, resulting in only two observable quartet peaks.")

    # Step 5: Final assignments and name
    print("\nStep 5: Final Assignment and Naming")
    print("The assignments for the proposed structure are:")
    print(f"{'145 ppm (s):':<15} Quaternary C of alkene (>C=)")
    print(f"{'112 ppm (t):':<15} Methylene C of alkene (=CH2)")
    print(f"{'48 ppm (t):':<15} Methylene C of the ethyl group (-CH2-)")
    print(f"{'27 ppm (d):':<15} Methine C (-CH-)")
    print(f"{'22 & 21 ppm (q):':<15} Three methyl C's (-CH3), with two overlapping")
    
    final_iupac_name = "2-ethyl-3-methyl-1-butene"
    print("\nTherefore, the IUPAC name of the compound is:")
    print(final_iupac_name)

if __name__ == '__main__':
    solve_nmr_puzzle()