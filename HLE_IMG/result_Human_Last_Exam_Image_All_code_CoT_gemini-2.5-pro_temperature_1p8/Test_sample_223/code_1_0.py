def calculate_molecular_formula():
    """
    Calculates the molecular formula of the organic cation of compound B based on the provided reaction scheme.
    """

    # Step 1: Define the molecular formula for the starting xanthenylium cation 'I' (C23H21O5+)
    cation_I = {'C': 23, 'H': 21, 'O': 5}
    print("Step 1: Determine the formula of the starting cation 'I' (1,8-dimethoxy-9-(2,6-dimethoxyphenyl)xanthen-9-ylium).")
    print("Formula of cation 'I' is C23H21O5(+).")
    print(f"  - Carbons: {cation_I['C']}")
    print(f"  - Hydrogens: {cation_I['H']}")
    print(f"  - Oxygens: {cation_I['O']}")
    print("-" * 30)

    # Step 2: Define the molecular formula for the reagent, methyl-3-aminopropionate (C4H9NO2)
    amine_reagent = {'C': 4, 'H': 9, 'N': 1, 'O': 2}
    print("Step 2: Determine the formula of the amine reagent (methyl-3-aminopropionate).")
    print("Formula of the amine is C4H9NO2.")
    print(f"  - Carbons: {amine_reagent['C']}")
    print(f"  - Hydrogens: {amine_reagent['H']}")
    print(f"  - Nitrogens: {amine_reagent['N']}")
    print(f"  - Oxygens: {amine_reagent['O']}")
    print("-" * 30)

    # Step 3: Identify the atoms lost during the reaction (H2O)
    lost_molecule = {'H': 2, 'O': 1}
    print("Step 3: Analyze the transformation.")
    print("The reaction forms an acridinium ring system with the elimination of one water molecule (H2O).")
    print(f"  - Hydrogens lost: {lost_molecule['H']}")
    print(f"  - Oxygens lost: {lost_molecule['O']}")
    print("-" * 30)

    # Step 4: Calculate the molecular formula of the final product cation B+
    print("Step 4: Calculate the molecular formula for the cation of compound B.")
    
    # Calculate Carbon atoms
    final_C = cation_I['C'] + amine_reagent['C']
    print(f"  - Final Carbons = (C from I) + (C from amine) = {cation_I['C']} + {amine_reagent['C']} = {final_C}")

    # Calculate Hydrogen atoms
    final_H = cation_I['H'] + amine_reagent['H'] - lost_molecule['H']
    print(f"  - Final Hydrogens = (H from I) + (H from amine) - (H lost) = {cation_I['H']} + {amine_reagent['H']} - {lost_molecule['H']} = {final_H}")
    
    # Calculate Nitrogen atoms
    final_N = amine_reagent['N']
    print(f"  - Final Nitrogens = (N from amine) = {amine_reagent['N']}")

    # Calculate Oxygen atoms
    final_O = cation_I['O'] + amine_reagent['O'] - lost_molecule['O']
    print(f"  - Final Oxygens = (O from I) + (O from amine) - (O lost) = {cation_I['O']} + {amine_reagent['O']} - {lost_molecule['O']} = {final_O}")
    print("-" * 30)

    # Step 5: Print the final molecular formula
    final_formula = f"C{final_C}H{final_H}NO{final_O}"
    print("The molecular formula of the organic cation of compound B is:")
    print(final_formula)

calculate_molecular_formula()