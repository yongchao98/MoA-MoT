def calculate_molecular_formula():
    """
    Calculates the molecular formula of compound B based on the reaction scheme.
    """
    # Step 1: Define the atomic composition of the starting cation 'I' (C23H21O5)
    cation_I = {'C': 23, 'H': 21, 'O': 5, 'N': 0}

    # Step 2: Define the atomic composition of the reagent, methyl-3-aminopropionate (C4H9NO2)
    reagent_B = {'C': 4, 'H': 9, 'N': 1, 'O': 2}

    # Step 3: Define the atomic composition of the byproduct, water (H2O)
    water = {'C': 0, 'H': 2, 'O': 1, 'N': 0}

    # Step 4: Calculate the composition of the product cation 'B'
    # B = I + reagent - water
    cation_B = {}
    all_elements = set(cation_I.keys()) | set(reagent_B.keys()) | set(water.keys())
    for element in all_elements:
        cation_B[element] = cation_I.get(element, 0) + reagent_B.get(element, 0) - water.get(element, 0)

    # Step 5: Format and print the final molecular formula
    c = cation_B.get('C', 0)
    h = cation_B.get('H', 0)
    n = cation_B.get('N', 0)
    o = cation_B.get('O', 0)

    print("The molecular formula of the organic cation of compound B is derived as follows:")
    print(f"Starting Cation 'I': C{cation_I['C']}H{cation_I['H']}O{cation_I['O']}")
    print(f"Reagent: C{reagent_B['C']}H{reagent_B['H']}N{reagent_B['N']}O{reagent_B['O']}")
    print(f"Byproduct: H{water['H']}O{water['O']}")
    print("\nCalculation:")
    print(f"Carbon atoms = {cation_I['C']} + {reagent_B['C']} - {water['C']} = {c}")
    print(f"Hydrogen atoms = {cation_I['H']} + {reagent_B['H']} - {water['H']} = {h}")
    print(f"Nitrogen atoms = {cation_I['N']} + {reagent_B['N']} - {water['N']} = {n}")
    print(f"Oxygen atoms = {cation_I['O']} + {reagent_B['O']} - {water['O']} = {o}")

    # Build the final formula string in standard order (C, H, then alphabetical)
    formula = f"C{c}H{h}"
    if n > 1:
        formula += f"N{n}"
    elif n == 1:
        formula += "N"
    if o > 1:
        formula += f"O{o}"
    elif o == 1:
        formula += "O"
        
    print(f"\nThe final molecular formula of the organic cation of compound B is: {formula}")

calculate_molecular_formula()