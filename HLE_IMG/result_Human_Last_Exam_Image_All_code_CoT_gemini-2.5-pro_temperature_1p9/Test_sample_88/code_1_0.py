import re

def parse_formula(formula):
    """Parses a chemical formula string into a dictionary of atom counts."""
    pattern = r'([A-Z][a-z]*)(\d*)'
    atoms = {}
    for element, count in re.findall(pattern, formula):
        atoms[element] = atoms.get(element, 0) + (int(count) if count else 1)
    return atoms

def add_formulas(dict1, dict2):
    """Adds the atom counts from two formula dictionaries."""
    result = dict1.copy()
    for element, count in dict2.items():
        result[element] = result.get(element, 0) + count
    return result

def format_formula_with_numbers(parsed_formula):
    """Formats a parsed formula dictionary into a string showing numbers, e.g., C9 H14 N2 O2"""
    parts = []
    # Ensure consistent order for printing
    for element in sorted(parsed_formula.keys()):
        count = parsed_formula[element]
        parts.append(f"{element}{count}")
    return " ".join(parts)


# --- Main script ---

# Molecular formulas from the problem analysis
sm_formula_str = "C9H14N2O2"
mp_formula_str = "C4H4O2"
product_b_formula_str = "C12H14N2O3"
methanol_formula_str = "CH4O"

# Parse formulas into dictionaries
sm_atoms = parse_formula(sm_formula_str)
mp_atoms = parse_formula(mp_formula_str)
product_b_atoms = parse_formula(product_b_formula_str)
methanol_atoms = parse_formula(methanol_formula_str)

# Calculate total atoms for reactants and products
reactant_atoms = add_formulas(sm_atoms, mp_atoms)
product_side_atoms = add_formulas(product_b_atoms, methanol_atoms)

# Print the analysis and verification
print("Analysis for the formation of Product B")
print("-" * 40)
print("Proposed Reaction Pathway:")
print(f"Starting Material ({sm_formula_str}) + Methyl Propiolate ({mp_formula_str}) -> Product B ({product_b_formula_str}) + Methanol ({methanol_formula_str})")
print("\nVerifying the atom balance for this reaction:")

# Print the equation with full formulas
print(f"\n{sm_formula_str} + {mp_formula_str} --> {product_b_formula_str} + {methanol_formula_str}\n")

print(f"Total atoms in Reactants: {sorted(reactant_atoms.items())}")
print(f"Total atoms in Products:  {sorted(product_side_atoms.items())}")

# Check if the reaction is balanced
if reactant_atoms == product_side_atoms:
    print("\nConclusion: The reaction is balanced. The proposed formation of Product B with elimination of methanol is stoichiometrically correct.")
else:
    print("\nConclusion: The reaction is NOT balanced.")

# Final requested output format
print("\nFinal equation with each number explicitly shown:")
print(f"1 * ({format_formula_with_numbers(sm_atoms)}) + 1 * ({format_formula_with_numbers(mp_atoms)}) -> 1 * ({format_formula_with_numbers(product_b_atoms)}) + 1 * ({format_formula_with_numbers(methanol_atoms)})")

# Calculate and print the molecular weight of Product B
atomic_weights = {'C': 12.011, 'H': 1.008, 'N': 14.007, 'O': 15.999}
mw_b = sum(atomic_weights[element] * count for element, count in product_b_atoms.items())
print(f"\nThe calculated molecular weight of Product B (C12H14N2O3) is: {mw_b:.2f} g/mol")