import collections

def molecular_calculation_and_verification():
    """
    Verifies the molecular formulas of products A, B, and C based on
    the hypothesized reaction pathways from the starting material.
    """

    # Define atomic weights for molecular weight calculation
    atomic_weights = {'C': 12.011, 'H': 1.008, 'N': 14.007, 'O': 15.999}

    # Define molecular formulas of all relevant species as dictionaries
    # Starting Material: (3,4-dihydro-2H-pyrrol-5-yl)proline
    sm = collections.Counter({'C': 9, 'H': 14, 'N': 2, 'O': 2})
    
    # Reagents and Fragments
    methyl_propiolate = collections.Counter({'C': 4, 'H': 4, 'O': 2})
    # Net change for mixed anhydride formation: + C2H2O
    anhydride_part = collections.Counter({'C': 2, 'H': 2, 'O': 1})
    methanol = collections.Counter({'C': 1, 'H': 4, 'O': 1})
    co2 = collections.Counter({'C': 1, 'O': 2})
    
    # Given product formulas
    product_a_formula = collections.Counter({'C': 14, 'H': 20, 'N': 2, 'O': 3})
    product_b_formula = collections.Counter({'C': 12, 'H': 14, 'N': 2, 'O': 3})
    product_c_formula = collections.Counter({'C': 11, 'H': 16, 'N': 2, 'O': 3})

    def get_formula_str(formula_dict):
        # Helper function to format a formula dictionary into a readable string
        return " ".join([f"{elem}{formula_dict[elem]}" for elem in ['C', 'H', 'N', 'O'] if elem in formula_dict and formula_dict[elem] > 0])
        
    def calculate_mw(formula_dict):
        # Helper function to calculate molecular weight
        mw = sum(atomic_weights[atom] * count for atom, count in formula_dict.items())
        return mw

    print("--- Molecular Formula Verification ---\n")
    
    # --- Verify Product C ---
    print(f"Verifying Product C ({get_formula_str(product_c_formula)})")
    # C is formed by creating a mixed anhydride from the starting material (SM)
    # SM (R-COOH) + Ac2O -> Product C (R-COO-Ac) + AcOH
    # The net change in the molecule is + C2H2O
    calculated_c = sm + anhydride_part
    print(f"Equation: SM + C2H2O = Product C")
    print(f"({get_formula_str(sm)}) + ({get_formula_str(anhydride_part)}) = {get_formula_str(calculated_c)}")
    if calculated_c == product_c_formula:
        print("Result: Calculation matches Product C formula.\n")
    else:
        print("Result: Mismatch.\n")
        
    # --- Verify Product B ---
    print(f"Verifying Product B ({get_formula_str(product_b_formula)})")
    # B is formed by addition of methyl propiolate (MP) to SM, then loss of methanol
    calculated_b = sm + methyl_propiolate - methanol
    print("Equation: SM + Methyl Propiolate - Methanol = Product B")
    print(f"({get_formula_str(sm)}) + ({get_formula_str(methyl_propiolate)}) - ({get_formula_str(methanol)}) = {get_formula_str(calculated_b)}")
    if calculated_b == product_b_formula:
        print("Result: Calculation matches Product B formula.\n")
    else:
        print("Result: Mismatch.\n")

    # --- Verify Product A ---
    print(f"Verifying Product A ({get_formula_str(product_a_formula)})")
    # A is formed by reaction of Product C with MP, followed by loss of CO2
    calculated_a = product_c_formula + methyl_propiolate - co2
    print("Equation: Product C + Methyl Propiolate - CO2 = Product A")
    print(f"({get_formula_str(product_c_formula)}) + ({get_formula_str(methyl_propiolate)}) - ({get_formula_str(co2)}) = {get_formula_str(calculated_a)}")
    if calculated_a == product_a_formula:
        print("Result: Calculation matches Product A formula.\n")
    else:
        print("Result: Mismatch.\n")

    # --- Final Molecular Weight Calculations ---
    print("--- Calculated Molecular Weights ---")
    print(f"Starting Material ({get_formula_str(sm)}): {calculate_mw(sm):.2f} g/mol")
    print(f"Product A ({get_formula_str(product_a_formula)}): {calculate_mw(product_a_formula):.2f} g/mol")
    print(f"Product B ({get_formula_str(product_b_formula)}): {calculate_mw(product_b_formula):.2f} g/mol")
    print(f"Product C ({get_formula_str(product_c_formula)}): {calculate_mw(product_c_formula):.2f} g/mol")

molecular_calculation_and_verification()