def calculate_molecular_formula_B():
    """
    Calculates the molecular formula of compound B based on the reaction scheme.
    """
    # Step 1: Define the atom counts for the known molecules based on the analysis.

    # Reagent for reaction B: methyl-3-aminopropionate (C4H9NO2)
    amine_B = {'C': 4, 'H': 9, 'N': 1, 'O': 2}

    # Starting xanthylium cation R (C25H20O5), as deduced from the reaction to form A
    cation_R = {'C': 25, 'H': 20, 'O': 5, 'N': 0}

    # Molecule eliminated during condensation: Water (H2O)
    water = {'H': 2, 'O': 1, 'C': 0, 'N': 0}

    # Step 2: Calculate the atom counts for the cation of compound B.
    # The reaction is: R+ + Amine -> B+ + H2O
    # So, B+ = R+ + Amine - H2O

    # Carbon count for B+
    c_b = cation_R['C'] + amine_B['C']
    # Hydrogen count for B+
    h_b = cation_R['H'] + amine_B['H'] - water['H']
    # Nitrogen count for B+
    n_b = cation_R['N'] + amine_B['N']
    # Oxygen count for B+
    o_b = cation_R['O'] + amine_B['O'] - water['O']

    # Step 3: Print the calculation steps and the final molecular formula.
    print("Calculating the molecular formula for the cation of compound B.")
    print("The reaction is a condensation of the starting cation R (C25H20O5+) with methyl-3-aminopropionate (C4H9NO2), eliminating one molecule of water (H2O).")
    print("\nCalculation steps for each element:")
    print(f"Carbon (C): {cation_R['C']} (from R) + {amine_B['C']} (from amine) = {c_b}")
    print(f"Hydrogen (H): {cation_R['H']} (from R) + {amine_B['H']} (from amine) - {water['H']} (from water) = {h_b}")
    print(f"Nitrogen (N): {cation_R['N']} (from R) + {amine_B['N']} (from amine) = {n_b}")
    print(f"Oxygen (O): {cation_R['O']} (from R) + {amine_B['O']} (from amine) - {water['O']} (from water) = {o_b}")

    molecular_formula_B = f"C{c_b}H{h_b}NO{o_b}"
    if n_b == 1:
        molecular_formula_B = f"C{c_b}H{h_b}NO{o_b}"
    else:
        molecular_formula_B = f"C{c_b}H{h_b}N{n_b}O{o_b}"

    print(f"\nThe molecular formula of the cation of compound B is: {molecular_formula_B}")

# Execute the function
calculate_molecular_formula_B()