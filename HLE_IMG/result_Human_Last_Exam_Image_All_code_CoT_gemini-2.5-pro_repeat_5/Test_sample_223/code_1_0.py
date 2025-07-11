def calculate_molecular_formula():
    """
    Calculates the molecular formula of compound B based on the reaction scheme.
    """
    # Step 1: Define the molecular formula of the starting xanthenium cation (C23H21O5+)
    start_cation = {'C': 23, 'H': 21, 'N': 0, 'O': 5}

    # Step 2: Define the molecular formula of the reagent, methyl-3-aminopropionate (C4H9NO2)
    reagent_amine = {'C': 4, 'H': 9, 'N': 1, 'O': 2}

    # Step 3: Define the molecular formula of the eliminated water molecule (H2O)
    eliminated_water = {'C': 0, 'H': 2, 'N': 0, 'O': 1}

    # Step 4: Calculate the formula for the cation of compound B
    product_cation_b = {}
    
    print("Calculating the molecular formula for the cation of Compound B:")
    
    # Carbon calculation
    c_final = start_cation['C'] + reagent_amine['C'] - eliminated_water['C']
    print(f"C: {start_cation['C']} + {reagent_amine['C']} - {eliminated_water['C']} = {c_final}")
    product_cation_b['C'] = c_final

    # Hydrogen calculation
    h_final = start_cation['H'] + reagent_amine['H'] - eliminated_water['H']
    print(f"H: {start_cation['H']} + {reagent_amine['H']} - {eliminated_water['H']} = {h_final}")
    product_cation_b['H'] = h_final

    # Nitrogen calculation
    n_final = start_cation['N'] + reagent_amine['N'] - eliminated_water['N']
    print(f"N: {start_cation['N']} + {reagent_amine['N']} - {eliminated_water['N']} = {n_final}")
    product_cation_b['N'] = n_final

    # Oxygen calculation
    o_final = start_cation['O'] + reagent_amine['O'] - eliminated_water['O']
    print(f"O: {start_cation['O']} + {reagent_amine['O']} - {eliminated_water['O']} = {o_final}")
    product_cation_b['O'] = o_final

    # Construct the final formula string
    # Handling the case where N is 1 (don't write N1)
    n_str = "N" if product_cation_b['N'] == 1 else f"N{product_cation_b['N']}"

    final_formula = f"C{product_cation_b['C']}H{product_cation_b['H']}{n_str}O{product_cation_b['O']}"
    
    print("\nThe molecular formula of the cation of compound B is:")
    print(final_formula)

calculate_molecular_formula()