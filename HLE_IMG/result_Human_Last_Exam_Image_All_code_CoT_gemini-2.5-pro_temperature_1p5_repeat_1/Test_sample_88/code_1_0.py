def print_formula(name, formula_dict):
    """Helper function to print a chemical formula nicely."""
    return f"{name} (C{formula_dict['C']}H{formula_dict['H']}N{formula_dict['N']}O{formula_dict['O']})"

def solve_chemistry_puzzle():
    """
    Solves the organic chemistry puzzle by calculating molecular formulas
    based on proposed reaction pathways.
    """
    # --- Define Known Molecular Formulas ---
    product_A_formula = {'C': 14, 'H': 20, 'N': 2, 'O': 3}
    product_B_formula = {'C': 12, 'H': 14, 'N': 2, 'O': 3}
    product_C_formula = {'C': 11, 'H': 16, 'N': 2, 'O': 3}

    # Reagents and Fragments
    methyl_propiolate = {'C': 4, 'H': 4, 'N': 0, 'O': 2}
    # Acetylation (from Ac2O) effectively replaces an H with a COCH3 group, net change is C2H2O
    acetyl_group_addition = {'C': 2, 'H': 2, 'N': 0, 'O': 1} 
    co2 = {'C': 1, 'H': 0, 'N': 0, 'O': 2}
    methyl_acetate = {'C': 3, 'H': 6, 'N': 0, 'O': 2}

    print("--- Analysis of the Reaction ---")
    
    # --- Step 1: Determine Starting Material (SM) from Product C ---
    print("\nStep 1: Determining the formula of the Starting Material (SM)")
    print(f"Hypothesis: Product C is the N-acetylated starting material.")
    print(f"Product C = SM + Acetyl_Group")
    print(f"Therefore, SM = Product C - Acetyl_Group")

    sm_formula = {
        'C': product_C_formula['C'] - acetyl_group_addition['C'],
        'H': product_C_formula['H'] - acetyl_group_addition['H'],
        'N': product_C_formula['N'] - acetyl_group_addition['N'],
        'O': product_C_formula['O'] - acetyl_group_addition['O'],
    }
    
    print("\nCalculation for Starting Material (SM):")
    print(f"  C: {product_C_formula['C']} - {acetyl_group_addition['C']} = {sm_formula['C']}")
    print(f"  H: {product_C_formula['H']} - {acetyl_group_addition['H']} = {sm_formula['H']}")
    print(f"  N: {product_C_formula['N']} - {acetyl_group_addition['N']} = {sm_formula['N']}")
    print(f"  O: {product_C_formula['O']} - {acetyl_group_addition['O']} = {sm_formula['O']}")
    print(f"Result: The starting material has the formula {print_formula('SM', sm_formula)}")

    # --- Step 2: Explain Formation of Product A ---
    print("\nStep 2: Proposing a pathway for Product A")
    print("Hypothesis: Product A is formed from N-acetyl SM (Product C) via decarboxylation to an azomethine ylide (AY), followed by cycloaddition with methyl propiolate (MP).")
    print("Pathway: C -> AY + CO2; AY + MP -> A")
    
    # Calculate formula of the azomethine ylide (AY)
    ay_formula = {
        'C': product_C_formula['C'] - co2['C'],
        'H': product_C_formula['H'] - co2['H'],
        'N': product_C_formula['N'] - co2['N'],
        'O': product_C_formula['O'] - co2['O'],
    }
    print(f"Intermediate {print_formula('AY', ay_formula)} is formed.")

    # Calculate formula of the final product A
    calc_A = {
        'C': ay_formula['C'] + methyl_propiolate['C'],
        'H': ay_formula['H'] + methyl_propiolate['H'],
        'N': ay_formula['N'] + methyl_propiolate['N'],
        'O': ay_formula['O'] + methyl_propiolate['O'],
    }

    print("\nCalculation for Product A:")
    print(f"Final Equation: {print_formula('AY', ay_formula)} + {print_formula('MP', methyl_propiolate)} -> Product A")
    print(f"  C: {ay_formula['C']} + {methyl_propiolate['C']} = {calc_A['C']}")
    print(f"  H: {ay_formula['H']} + {methyl_propiolate['H']} = {calc_A['H']}")
    print(f"  N: {ay_formula['N']} + {methyl_propiolate['N']} = {calc_A['N']}")
    print(f"  O: {ay_formula['O']} + {methyl_propiolate['O']} = {calc_A['O']}")
    print(f"Result: Calculated formula is {print_formula('', calc_A)}, which matches Product A.")

    # --- Step 3: Explain Formation of Product B ---
    print("\nStep 3: Proposing a pathway for Product B")
    print("Hypothesis: Product B is formed from N-acetyl SM (Product C) reacting with methyl propiolate (MP), followed by elimination of methyl acetate (MeOAc).")
    print("Pathway: C + MP -> [Intermediate] -> B + MeOAc")

    calc_B = {
        'C': product_C_formula['C'] + methyl_propiolate['C'] - methyl_acetate['C'],
        'H': product_C_formula['H'] + methyl_propiolate['H'] - methyl_acetate['H'],
        'N': product_C_formula['N'] + methyl_propiolate['N'] - methyl_acetate['N'],
        'O': product_C_formula['O'] + methyl_propiolate['O'] - methyl_acetate['O'],
    }

    print("\nCalculation for Product B:")
    print(f"Final Equation: {print_formula('Product C', product_C_formula)} + {print_formula('MP', methyl_propiolate)} - {print_formula('MeOAc', methyl_acetate)} -> Product B")
    print(f"  C: {product_C_formula['C']} + {methyl_propiolate['C']} - {methyl_acetate['C']} = {calc_B['C']}")
    print(f"  H: {product_C_formula['H']} + {methyl_propiolate['H']} - {methyl_acetate['H']} = {calc_B['H']}")
    print(f"  N: {product_C_formula['N']} + {methyl_propiolate['N']} - {methyl_acetate['N']} = {calc_B['N']}")
    print(f"  O: {product_C_formula['O']} + {methyl_propiolate['O']} - {methyl_acetate['O']} = {calc_B['O']}")
    print(f"Result: Calculated formula is {print_formula('', calc_B)}, which matches Product B.")

solve_chemistry_puzzle()
>>> 