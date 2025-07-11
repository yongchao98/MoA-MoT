import re

def calculate_molar_mass(formula):
    """
    Calculates the molar mass of a chemical compound given its molecular formula.
    Also returns a string representing the calculation steps.
    """
    atomic_weights = {
        'C': 12.011,
        'H': 1.008,
        'N': 14.007,
        'O': 15.999
    }

    # Regex to parse the formula, e.g., C14H20N2O3
    # It finds an uppercase letter (element) followed by optional digits (count)
    atom_counts = re.findall(r'([A-Z])(\d*)', formula)

    total_mass = 0
    calculation_str_parts = []
    
    for element, count_str in atom_counts:
        count = int(count_str) if count_str else 1
        if element in atomic_weights:
            mass = count * atomic_weights[element]
            total_mass += mass
            calculation_str_parts.append(f"{count} * {atomic_weights[element]}")
        else:
            print(f"Warning: Atomic weight for element {element} not found.")
            return None, None

    calculation_str = " + ".join(calculation_str_parts)
    return total_mass, calculation_str

# Molecular formulas of the products
products = {
    'A': 'C14H20N2O3',
    'B': 'C12H14N2O3',
    'C': 'C11H16N2O3'
}

# Calculate and print the molar mass for each product
print("Calculating the molar masses of the products:\n")

for product_name, formula in products.items():
    molar_mass, calc_str = calculate_molar_mass(formula)
    if molar_mass is not None:
        print(f"Product {product_name}:")
        print(f"  Formula: {formula}")
        print(f"  Calculation: M = {calc_str}")
        print(f"  Molar Mass = {molar_mass:.3f} g/mol\n")

# The major product in this reaction is A. The final answer will be its molar mass.
# Final Answer (Molar Mass of Product A): 264.325