def calculate_formula(base_formula, groups):
    """Calculates the final chemical formula from a base and attached groups."""
    final_formula = base_formula.copy()
    for group_formula, count in groups:
        for element, num in group_formula.items():
            final_formula[element] = final_formula.get(element, 0) + num * count
    return final_formula

def format_formula(formula_dict):
    """Formats a dictionary into a chemical formula string."""
    # Standard order for organic compounds: C, H, then alphabetical for others
    order = ['C', 'H']
    sorted_elements = sorted([el for el in formula_dict if el not in order])
    
    formula_str = ""
    for el in order:
        if el in formula_dict and formula_dict[el] > 0:
            formula_str += f"{el}{formula_dict[el] if formula_dict[el] > 1 else ''}"
    for el in sorted_elements:
        if el in formula_dict and formula_dict[el] > 0:
            formula_str += f"{el}{formula_dict[el] if formula_dict[el] > 1 else ''}"
    return formula_str

# Define elemental compositions
phenyl_group_core = {'C': 6, 'H': 3} # C6H3 core for a 1,2,3-trisubstituted ring
methoxy_group = {'C': 1, 'H': 3, 'O': 1}
hydroxyl_group = {'O': 1, 'H': 1}
central_carbon_ion = {'C': 1}

# --- Starting Material: tris(2,3-dimethoxyphenyl)methylium ion ---
# Each arm is a C6H3 ring with two methoxy groups
dimethoxyphenyl_groups = [
    (phenyl_group_core, 3),
    (methoxy_group, 6)
]
# Calculate the formula for the three arms together
arms_formula = calculate_formula({}, dimethoxyphenyl_groups)

# Add the central carbon to get the final cation formula
start_material_formula = calculate_formula(central_carbon_ion, [(arms_formula, 1)])
start_material_name = "tris(2,3-dimethoxyphenyl)methylium ion"
start_formula_str = format_formula(start_material_formula)

print(f"Starting Material: {start_material_name}")
print(f"Formula: [{start_formula_str}]+")
print("-" * 30)

# --- Product A: tris(2,3-dihydroxyphenyl)methylium ion ---
# Reaction: 6 x OCH3 are hydrolyzed to 6 x OH
dihydroxyphenyl_groups = [
    (phenyl_group_core, 3),
    (hydroxyl_group, 6)
]
# Calculate the formula for the three arms of the product
product_arms_formula = calculate_formula({}, dihydroxyphenyl_groups)

# Add the central carbon to get the final cation formula
product_formula = calculate_formula(central_carbon_ion, [(product_arms_formula, 1)])
product_name = "tris(2,3-dihydroxyphenyl)methylium ion"
product_formula_str = format_formula(product_formula)

print(f"Compound A: {product_name}")
print(f"Formula: [{product_formula_str}]+")

# The reaction is: [C25H27O6]+ + 6 H2O -> [C19H15O6]+ + 6 CH3OH
print("\nOverall Reaction:")
print(f"[{start_formula_str}]+ + 6H2O --(0.1 M HCl, reflux)--> [{product_formula_str}]+ + 6CH3OH")
