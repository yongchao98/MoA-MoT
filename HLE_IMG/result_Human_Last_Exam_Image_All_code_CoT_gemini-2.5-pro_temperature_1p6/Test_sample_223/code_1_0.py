def calculate_molecular_formula():
    """
    This function calculates the molecular formula of compound B based on the reaction provided.
    The reaction involves a xanthylium cation reacting with methyl-3-aminopropionate to form
    an acridinium cation (B) and water.
    """
    
    # Molecular formula of the starting xanthylium cation: [C23H21O5]+
    xanthylium_cation = {'C': 23, 'H': 21, 'O': 5, 'N': 0}

    # Molecular formula of methyl-3-aminopropionate: C4H9NO2
    amine = {'C': 4, 'H': 9, 'N': 1, 'O': 2}

    # Molecular formula of water, the byproduct: H2O
    water = {'C': 0, 'H': 2, 'O': 1, 'N': 0}

    # Calculate the molecular formula of the product cation B
    # Formula(B+) = Formula(xanthylium+) + Formula(amine) - Formula(water)
    b_cation = {}
    
    # Calculate atom counts
    b_cation['C'] = xanthylium_cation['C'] + amine['C'] - water['C']
    b_cation['H'] = xanthylium_cation['H'] + amine['H'] - water['H']
    b_cation['N'] = xanthylium_cation['N'] + amine['N'] - water['N']
    b_cation['O'] = xanthylium_cation['O'] + amine['O'] - water['O']

    # Print the explanation and the calculation steps
    print("To find the molecular formula for the cation of compound B, we use the following reaction stoichiometry:")
    print("Formula(B cation) = Formula(Starting Cation) + Formula(Amine) - Formula(Water)\n")
    
    print("Starting Cation Formula: C23 H21 O5")
    print("Amine Formula: C4 H9 N1 O2")
    print("Water Formula: H2 O1\n")

    print("Calculation for each element:")
    print(f"Carbon (C): {xanthylium_cation['C']} + {amine['C']} = {b_cation['C']}")
    print(f"Hydrogen (H): {xanthylium_cation['H']} + {amine['H']} - {water['H']} = {b_cation['H']}")
    print(f"Nitrogen (N): {xanthylium_cation['N']} + {amine['N']} - {water['N']} = {b_cation['N']}")
    print(f"Oxygen (O): {xanthylium_cation['O']} + {amine['O']} - {water['O']} = {b_cation['O']}\n")
    
    # Construct the final formula string
    final_formula = f"C{b_cation['C']}H{b_cation['H']}N{b_cation['N']}O{b_cation['O']}"

    print(f"The molecular formula of the cation in compound B is: {final_formula}")

# Run the calculation
calculate_molecular_formula()