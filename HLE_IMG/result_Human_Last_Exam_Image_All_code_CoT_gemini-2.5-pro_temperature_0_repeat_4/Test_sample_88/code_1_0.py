def calculate_formula():
    """
    Calculates the molecular formulas of products A, B, and C based on the reaction scheme.
    """
    # Define molecular formulas as dictionaries
    sm = {'C': 9, 'H': 14, 'N': 2, 'O': 2} # Starting Material: (3,4-dihydro-2H-pyrrol-5-yl)proline
    acetyl_addition = {'C': 2, 'H': 2, 'O': 1} # Net change for acetylation (-H, +COCH3)
    methyl_propiolate = {'C': 4, 'H': 4, 'O': 2}
    co2 = {'C': 1, 'H': 0, 'N': 0, 'O': 2}
    c2h6_fragment = {'C': 2, 'H': 6, 'N': 0, 'O': 0}

    # Helper function to format formula string
    def format_formula(formula_dict, name):
        return f"Product {name}: C{formula_dict['C']}H{formula_dict['H']}N{formula_dict['N']}O{formula_dict['O']}"

    # --- Calculate Product C ---
    print("Calculating the molecular formula for Product C...")
    print(f"The starting material (SM) has the formula: C{sm['C']}H{sm['H']}N{sm['N']}O{sm['O']}")
    print("Product C is formed by N-acetylation of the starting material.")
    print(f"This corresponds to a net addition of a C{acetyl_addition['C']}H{acetyl_addition['H']}O{acetyl_addition['O']} group.")
    
    product_c = {
        'C': sm['C'] + acetyl_addition['C'],
        'H': sm['H'] + acetyl_addition['H'],
        'N': sm['N'],
        'O': sm['O'] + acetyl_addition['O']
    }
    
    print(f"Calculation: C = {sm['C']} + {acetyl_addition['C']} = {product_c['C']}")
    print(f"             H = {sm['H']} + {acetyl_addition['H']} = {product_c['H']}")
    print(f"             N = {sm['N']}")
    print(f"             O = {sm['O']} + {acetyl_addition['O']} = {product_c['O']}")
    print(f"The calculated formula for Product C is C{product_c['C']}H{product_c['H']}N{product_c['N']}O{product_c['O']}.\n")

    # --- Calculate Product A ---
    print("Calculating the molecular formula for Product A...")
    print(f"Product A is formed from Product C (C{product_c['C']}H{product_c['H']}N{product_c['N']}O{product_c['O']}) reacting with methyl propiolate (C{methyl_propiolate['C']}H{methyl_propiolate['H']}O{methyl_propiolate['O']}) and losing CO2 (C{co2['C']}O{co2['O']}).")
    
    product_a = {
        'C': product_c['C'] + methyl_propiolate['C'] - co2['C'],
        'H': product_c['H'] + methyl_propiolate['H'],
        'N': product_c['N'],
        'O': product_c['O'] + methyl_propiolate['O'] - co2['O']
    }

    print(f"Calculation: C = {product_c['C']} + {methyl_propiolate['C']} - {co2['C']} = {product_a['C']}")
    print(f"             H = {product_c['H']} + {methyl_propiolate['H']} = {product_a['H']}")
    print(f"             N = {product_c['N']}")
    print(f"             O = {product_c['O']} + {methyl_propiolate['O']} - {co2['O']} = {product_a['O']}")
    print(f"The calculated formula for Product A is C{product_a['C']}H{product_a['H']}N{product_a['N']}O{product_a['O']}.\n")

    # --- Calculate Product B ---
    print("Calculating the molecular formula for Product B...")
    print(f"Product B can be derived from Product A (C{product_a['C']}H{product_a['H']}N{product_a['N']}O{product_a['O']}) by the loss of a C{c2h6_fragment['C']}H{c2h6_fragment['H']} fragment.")

    product_b = {
        'C': product_a['C'] - c2h6_fragment['C'],
        'H': product_a['H'] - c2h6_fragment['H'],
        'N': product_a['N'],
        'O': product_a['O']
    }
    
    print(f"Calculation: C = {product_a['C']} - {c2h6_fragment['C']} = {product_b['C']}")
    print(f"             H = {product_a['H']} - {c2h6_fragment['H']} = {product_b['H']}")
    print(f"             N = {product_a['N']}")
    print(f"             O = {product_a['O']}")
    print(f"The calculated formula for Product B is C{product_b['C']}H{product_b['H']}N{product_b['N']}O{product_b['O']}.\n")

    print("--- Summary of Final Formulas ---")
    print(format_formula(product_a, 'A'))
    print(format_formula(product_b, 'B'))
    print(format_formula(product_c, 'C'))

calculate_formula()