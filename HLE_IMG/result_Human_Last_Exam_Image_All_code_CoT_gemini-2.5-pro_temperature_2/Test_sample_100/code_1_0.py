def calculate_product_formula():
    """
    Calculates the molecular formula of the product through a three-step reaction.
    """
    # Step 0: Define the molecular formula of the starting material
    # C15H14F3NO2
    formula = {'C': 15, 'H': 14, 'F': 3, 'N': 1, 'O': 2}
    
    print("This script calculates the molecular formula of the final product in the reaction sequence.")
    print("-" * 30)
    
    def format_formula(f):
        # Format the dictionary into a string, e.g., C15H14F3NO2
        # Order: C, H, then alphabetically for the rest
        parts = []
        if 'C' in f and f['C'] > 0:
            parts.append(f"C{f['C'] if f['C'] > 1 else ''}")
        if 'H' in f and f['H'] > 0:
            parts.append(f"H{f['H'] if f['H'] > 1 else ''}")
        
        # Sort remaining keys alphabetically
        other_keys = sorted([k for k in f.keys() if k not in ['C', 'H']])
        
        for key in other_keys:
            if f[key] > 0:
                parts.append(f"{key}{f[key] if f[key] > 1 else ''}")
        return "".join(parts)

    print(f"Step 0: Starting Material Formula: {format_formula(formula)}")
    print("\n")

    # Store initial counts for the final equation
    initial_counts = formula.copy()
    
    # Step 1: Deprotection with CAN
    # Net change is removal of C8H8O
    change1 = {'C': -8, 'H': -8, 'O': -1}
    formula['C'] += change1.get('C', 0)
    formula['H'] += change1.get('H', 0)
    formula['O'] += change1.get('O', 0)
    print(f"Step 1: After deprotection (removal of C8H8O).")
    print(f"         Intermediate 1 Formula: {format_formula(formula)}")
    print("\n")

    # Step 2: Hydrogenation with Pd/C, H2
    # Net change is addition of H2
    change2 = {'H': 2}
    formula['H'] += change2.get('H', 0)
    print(f"Step 2: After hydrogenation (addition of H2).")
    print(f"         Intermediate 2 Formula: {format_formula(formula)}")
    print("\n")
    
    # Step 3: Hydrolysis with HCl
    # Net change is addition of H2O
    change3 = {'H': 2, 'O': 1}
    formula['H'] += change3.get('H', 0)
    formula['O'] += change3.get('O', 0)
    print(f"Step 3: After lactam hydrolysis (addition of H2O).")
    print(f"         Final Product Formula: {format_formula(formula)}")
    print("\n")

    # Final calculation breakdown as requested
    print("-" * 30)
    print("Final Atom Count Calculation:")
    
    final_C = initial_counts['C'] + change1.get('C', 0)
    print(f"Carbon (C):   {initial_counts['C']} {change1.get('C', 0):+d} = {final_C}")
    
    final_H = initial_counts['H'] + change1.get('H', 0) + change2.get('H', 0) + change3.get('H', 0)
    print(f"Hydrogen (H): {initial_counts['H']} {change1.get('H', 0):+d} {change2.get('H', 0):+d} {change3.get('H', 0):+d} = {final_H}")
    
    final_F = initial_counts['F']
    print(f"Fluorine (F): {initial_counts['F']} = {final_F}")

    final_N = initial_counts['N']
    print(f"Nitrogen (N): {initial_counts['N']} = {final_N}")
    
    final_O = initial_counts['O'] + change1.get('O', 0) + change3.get('O', 0)
    print(f"Oxygen (O):   {initial_counts['O']} {change1.get('O', 0):+d} {change3.get('O', 0):+d} = {final_O}")

    print("-" * 30)
    print(f"The final molecular formula of the product is {format_formula(formula)}.")

calculate_product_formula()
<<<C7H10F3NO2>>>