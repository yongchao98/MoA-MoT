def solve_chemistry_problem():
    """
    This script explains the reaction steps and calculates the molecular formula for compound C.
    """

    # Helper function to pretty-print a molecular formula from a dictionary
    def format_formula(formula_dict):
        # Standard order for organic compounds: C, H, then alphabetical
        order = ['C', 'H']
        atoms = sorted([atom for atom in formula_dict if formula_dict.get(atom, 0) > 0])
        
        formula_str = ""
        for atom in order:
            if atom in atoms:
                count = formula_dict[atom]
                formula_str += atom
                if count > 1:
                    formula_str += str(count)
                atoms.remove(atom)
        
        for atom in atoms:
            count = formula_dict[atom]
            formula_str += atom
            if count > 1:
                formula_str += str(count)
                
        return formula_str

    # Helper function to add two molecular formulas
    def add_formulas(f1, f2):
        res = f1.copy()
        for atom, count in f2.items():
            res[atom] = res.get(atom, 0) + count
        return res

    # Helper function to subtract one molecular formula from another
    def sub_formulas(f1, f2):
        res = f1.copy()
        for atom, count in f2.items():
            res[atom] = res.get(atom, 0) - count
        return res

    print("Step-by-step deduction of Compound C:")
    print("-" * 40)
    
    # Define molecular formulas of key fragments
    xanthylium_core = {'C': 13, 'H': 9, 'O': 1}
    methoxy_group = {'C': 1, 'H': 3, 'O': 1}
    diethylamino_group = {'C': 4, 'H': 10, 'N': 1}
    hydroxyl_group = {'O': 1, 'H': 1}

    print("1. Identifying the core structure:")
    print("The reaction creates a xanthylium cation core. Its molecular formula is C13H9O+.")
    print(f"   - Core composition: {format_formula(xanthylium_core)}")
    print("-" * 40)

    print("2. Determining Compound B (3,6-bis(diethylamino)-1,8-dimethoxyxanthylium):")
    print("   - Start with the xanthylium core (C13H9O+).")
    print("   - Replace 4 H atoms with 2 diethylamino groups and 2 methoxy groups.")
    
    # Calculation for B
    formula_b = sub_formulas(xanthylium_core, {'H': 4})
    formula_b = add_formulas(formula_b, diethylamino_group)
    formula_b = add_formulas(formula_b, diethylamino_group)
    formula_b = add_formulas(formula_b, methoxy_group)
    formula_b = add_formulas(formula_b, methoxy_group)
    print(f"   - Molecular formula of B is {format_formula(formula_b)}+.")
    print("-" * 40)
    
    print("3. Determining Compound C (3,6-bis(diethylamino)-1,8-dihydroxyxanthylium):")
    print("   - Start with Compound B.")
    print("   - The reaction with LiI replaces two methoxy groups (-OCH3) with two hydroxyl groups (-OH).")
    
    # Calculation for C from B
    formula_c = formula_b.copy()
    formula_c = sub_formulas(formula_c, methoxy_group)
    formula_c = sub_formulas(formula_c, methoxy_group)
    formula_c = add_formulas(formula_c, hydroxyl_group)
    formula_c = add_formulas(formula_c, hydroxyl_group)

    print(f"\nFinal Result for Compound C:")
    print("=" * 40)
    print(f"The name of Compound C is: 3,6-bis(diethylamino)-1,8-dihydroxyxanthylium cation")
    print(f"The molecular formula of the cation is: {format_formula(formula_c)}+")
    print("\nThe final 'equation' is the molecular formula. The numbers in this equation are the atom counts:")
    for atom, count in sorted(formula_c.items()):
        if count > 0:
            print(f"  Number of {atom} atoms: {count}")
    print("=" * 40)


if __name__ == "__main__":
    solve_chemistry_problem()
