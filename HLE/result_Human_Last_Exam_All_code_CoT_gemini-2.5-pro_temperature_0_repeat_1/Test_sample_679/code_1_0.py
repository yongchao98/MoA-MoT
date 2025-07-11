def solve_hydrocarbon_structure():
    """
    Determines the IUPAC name of a hydrocarbon from its molecular formula and 13C NMR data,
    and prints the step-by-step reasoning.
    """
    
    # Step 1: Analyze Molecular Formula
    print("Step 1: Analyzing the Molecular Formula")
    print("The molecular formula is C7H14.")
    print("For a saturated acyclic alkane with 7 carbons, the formula would be CnH(2n+2) = C7H16.")
    print("The given formula has two fewer hydrogens, which indicates one degree of unsaturation. This is likely a double bond or a ring.\n")

    # Step 2: Analyze 13C NMR Data
    print("Step 2: Analyzing the 13C NMR Data")
    print("The NMR data is: 145(s), 112(t), 48(t), 27(d), 22(q), 21(q).")
    print("Multiplicity key: (s) = quaternary C (0 H), (d) = CH (1 H), (t) = CH2 (2 H), (q) = CH3 (3 H).")
    print("- The signals at 145 ppm and 112 ppm are in the alkene region (typically 100-150 ppm), confirming a C=C double bond.")
    print("- The signal at 145(s) corresponds to a quaternary alkene carbon (>C=) with no attached hydrogens.")
    print("- The signal at 112(t) corresponds to a terminal alkene carbon (=CH2) with two attached hydrogens.")
    print("- These two signals together strongly indicate a >C=CH2 fragment in the molecule.\n")

    # Step 3: Reconcile Carbon and Hydrogen Counts
    print("Step 3: Reconciling the Data with the Formula")
    print("The formula has 7 carbons, but there are only 6 NMR signals. This means two carbons in the molecule are chemically equivalent and produce a single signal.")
    print("Let's count the hydrogens based on the multiplicities of the 6 signals:")
    print("H count = 0 (from 145s) + 2 (from 112t) + 2 (from 48t) + 1 (from 27d) + 3 (from 22q) + 3 (from 21q) = 11 hydrogens.")
    print("The molecular formula requires 14 hydrogens. The difference is 14 - 11 = 3 hydrogens.")
    print("The only way to account for both the missing 1 carbon and 3 hydrogens is if one of the quartet (q) signals represents two equivalent methyl (CH3) groups.\n")

    # Step 4: Assemble the Structure
    print("Step 4: Assembling the Molecular Structure")
    print("The structural fragments of the molecule are:")
    print("- A >C=CH2 group (from 145(s) and 112(t))")
    print("- A -CH2- group (from 48(t))")
    print("- A -CH- group (from 27(d))")
    print("- One non-equivalent -CH3 group (from one 'q' signal)")
    print("- Two equivalent -CH3 groups (from the other 'q' signal)")
    print("\nAssembly process:")
    print("- The two equivalent -CH3 groups must be attached to the -CH- group, forming an isopropyl group: -CH(CH3)2.")
    print("- The -CH2- group and the single -CH3 group combine to form an ethyl group: -CH2CH3.")
    print("- The quaternary carbon of the >C=CH2 unit has two available bonds. We attach the ethyl and isopropyl groups to it.")
    print("The final structure is: (CH3CH2)-C(=CH2)-CH(CH3)2\n")

    # Step 5: Determine the IUPAC Name
    print("Step 5: Determining the IUPAC Name")
    print("For the structure CH3-CH(CH3)-C(=CH2)-CH2-CH3:")
    print("1. The longest continuous carbon chain is a 5-carbon chain (pentane).")
    print("2. The double bond is outside the main chain (exocyclic), so it is named as a 'methylidene' substituent.")
    print("3. We number the chain to give the substituents the lowest possible numbers.")
    print("   Numbering from the end closer to the isopropyl fragment gives locants 2 (for the methyl group) and 3 (for the methylidene group).")
    print("   CH3(5)-CH2(4)-C(3)(=CH2)-CH(2)(CH3)-CH3(1)")
    print("4. The substituents are listed alphabetically (methyl, methylidene).")
    
    final_name = "2-methyl-3-methylidenepentane"
    print("\nThe IUPAC name of the compound is:")
    print(final_name)

solve_hydrocarbon_structure()