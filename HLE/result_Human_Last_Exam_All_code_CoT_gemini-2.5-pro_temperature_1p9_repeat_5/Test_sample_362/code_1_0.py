def get_formatted_formula(formula_dict):
    """
    Formats a dictionary of atom counts into a standard chemical formula string.
    Order: C, H, then other elements alphabetically.
    """
    # Define the standard order for organic compounds
    order = ['C', 'H']
    formula_str = ""
    
    # Append atoms in the defined order (C, H)
    for atom in order:
        if atom in formula_dict:
            count = formula_dict[atom]
            formula_str += atom
            if count > 1:
                formula_str += str(count)
    
    # Append remaining atoms in alphabetical order
    other_atoms = sorted([key for key in formula_dict if key not in order])
    for atom in other_atoms:
        count = formula_dict[atom]
        formula_str += atom
        if count > 1:
            formula_str += str(count)
            
    return formula_str

def main():
    """
    This script outlines the Wittig reaction of pivalaldehyde with a specified ylide,
    presenting the reactants, products, and the final balanced equation.
    """
    # Reactant 1: Pivalaldehyde
    reactant_1 = {
        "name": "Pivalaldehyde",
        "formula_dict": {'C': 5, 'H': 10, 'O': 1}
    }
    
    # Reactant 2: (2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane (the ylide)
    reactant_2 = {
        "name": "(2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane",
        "formula_dict": {'C': 26, 'H': 22, 'Cl': 1, 'P': 1}
    }

    # Product 1: The alkene
    product_1 = {
        "name": "1-(2-chlorophenyl)-4,4-dimethylpent-2-ene",
        "formula_dict": {'C': 13, 'H': 17, 'Cl': 1}
    }

    # Product 2: The byproduct
    product_2 = {
        "name": "Triphenylphosphine oxide",
        "formula_dict": {'C': 18, 'H': 15, 'O': 1, 'P': 1}
    }

    print("--- Wittig Reaction Analysis ---")
    print(f"\nReactant 1: {reactant_1['name']}")
    print(f"   Formula: {get_formatted_formula(reactant_1['formula_dict'])}")
    
    print(f"\nReactant 2: {reactant_2['name']}")
    print(f"   Formula: {get_formatted_formula(reactant_2['formula_dict'])}")

    print("\n      ||")
    print("      \\/")

    print(f"\nProduct 1 (Alkene): {product_1['name']}")
    print(f"   Formula: {get_formatted_formula(product_1['formula_dict'])}")

    print(f"\nProduct 2 (Byproduct): {product_2['name']}")
    print(f"   Formula: {get_formatted_formula(product_2['formula_dict'])}")
    
    print("\n--- Final Balanced Equation ---")
    # The instruction "output each number in the final equation" is interpreted as showing
    # the stoichiometric coefficient and the atom counts within each formula.
    r1_formula = get_formatted_formula(reactant_1['formula_dict'])
    r2_formula = get_formatted_formula(reactant_2['formula_dict'])
    p1_formula = get_formatted_formula(product_1['formula_dict'])
    p2_formula = get_formatted_formula(product_2['formula_dict'])

    print(f"1 {r1_formula} + 1 {r2_formula} -> 1 {p1_formula} + 1 {p2_formula}")

if __name__ == "__main__":
    main()