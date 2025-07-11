import re
from collections import defaultdict

def parse_formula(formula):
    """
    Parses a chemical formula string and returns a dictionary of its atom counts.
    For example, 'C9H14N2O2' becomes {'C': 9, 'H': 14, 'N': 2, 'O': 2}.
    """
    # This regex finds elements (e.g., C, Cl) followed by an optional number.
    pairs = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    atom_counts = defaultdict(int)
    for element, count_str in pairs:
        # If no number follows an element, its count is 1.
        count = int(count_str) if count_str else 1
        atom_counts[element] += count
    return atom_counts

def calculate_total_atoms(formulas):
    """
    Calculates the total atom counts for a list of chemical formulas.
    """
    total_counts = defaultdict(int)
    for formula in formulas:
        counts = parse_formula(formula)
        for element, count in counts.items():
            total_counts[element] += count
    return total_counts

def get_numbers_from_formulas(formulas):
    """
    Extracts all numerical values from a list of chemical formulas.
    If an element has no number, it is treated as 1.
    """
    all_numbers = []
    for formula in formulas:
        pairs = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
        for _, count_str in pairs:
            count = int(count_str) if count_str else 1
            all_numbers.append(count)
    return all_numbers

# Define the molecular formulas for the reactants and products in the proposed reaction for Product B.
# Starting Material: C9H14N2O2
# Methyl Propiolate: C4H4O2
# Product B: C12H14N2O3
# Methanol (eliminated molecule): CH4O
reactants_formulas = ["C9H14N2O2", "C4H4O2"]
products_formulas = ["C12H14N2O3", "CH4O"]

# Create the full chemical equation string for display.
equation_str = f"{' + '.join(reactants_formulas)} -> {' + '.join(products_formulas)}"
print(f"Proposed reaction for the formation of Product B:\n{equation_str}\n")

# Extract and print all numbers from the formulas in the equation.
all_formulas = reactants_formulas + products_formulas
numbers_in_equation = get_numbers_from_formulas(all_formulas)
print(f"The numbers in the final equation are: {', '.join(map(str, numbers_in_equation))}\n")

# Calculate the total atoms for the reactants and products side.
reactant_atoms = calculate_total_atoms(reactants_formulas)
product_atoms = calculate_total_atoms(products_formulas)

# Verify and print the atom balance for each element.
print("Verifying atom balance:")
all_elements = sorted(set(reactant_atoms.keys()) | set(product_atoms.keys()))

is_balanced = True
for element in all_elements:
    r_count = reactant_atoms.get(element, 0)
    p_count = product_atoms.get(element, 0)
    print(f"Element {element}: Reactants = {r_count}, Products = {p_count}")
    if r_count != p_count:
        is_balanced = False

# Print the final conclusion on whether the equation is balanced.
print("\n" + "="*25)
if is_balanced:
    print("Conclusion: The chemical equation is balanced.")
else:
    print("Conclusion: The chemical equation is NOT balanced.")
print("="*25)
