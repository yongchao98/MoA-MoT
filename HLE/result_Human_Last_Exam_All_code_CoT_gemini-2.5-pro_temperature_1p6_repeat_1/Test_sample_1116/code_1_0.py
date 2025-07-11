import collections

def calculate_molecular_weight(formula_dict, atomic_weights_dict):
    """Calculates the molecular weight from a dictionary of element counts."""
    total_weight = 0
    calculation_parts = []
    # Sort for consistent order
    for element in sorted(formula_dict.keys()):
        count = formula_dict[element]
        weight = atomic_weights_dict[element]
        total_weight += count * weight
        calculation_parts.append(f"{count} * {weight} ({element})")
    
    calculation_equation = " + ".join(calculation_parts)
    return total_weight, calculation_equation

# --- Analysis and Product Identification ---

# Atomic weights of relevant elements (g/mol)
ATOMIC_WEIGHTS = {
    'C': 12.011, 'H': 1.008, 'N': 14.007, 'O': 15.999, 'S': 32.06, 'Br': 79.904,
}

# The reaction replaces 2 Hydrogen atoms with 2 Bromine atoms.
# Molecular formula of the starting material: C35 H43 N O2 S4
# We derive the product's formula by subtracting 2 H and adding 2 Br.
product_formula = collections.OrderedDict([
    ('C', 35),
    ('H', 41),
    ('Br', 2),
    ('N', 1),
    ('O', 2),
    ('S', 4)
])

# Construct the product name based on the reaction
side_chain_brominated = "5-bromo-4-(2-ethylhexyl)thiophen-2-yl"
core = "5-methyl-4H-dithieno[3,2-e:2',3'-g]isoindole-4,6(5H)-dione"
product_name = f"2,8-bis({side_chain_brominated})-{core}"

# Calculate Molecular Weight
mw, mw_equation = calculate_molecular_weight(product_formula, ATOMIC_WEIGHTS)

# --- Output the Results ---
print("Based on the reaction analysis, the new spot is the dibrominated product.")
print("-" * 50)
print(f"Product Name: {product_name}")
print(f"\nMolecular Formula: C35H41Br2NO2S4")
print(f"\nMolecular Weight Calculation:")
print(f"Equation: {mw_equation}")
print(f"\nFinal Molecular Weight: {mw:.3f} g/mol")
print("-" * 50)
