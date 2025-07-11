import re
from collections import Counter

def parse_formula(formula: str) -> Counter:
    """Parses a molecular formula string into a Counter of its atoms."""
    atoms = Counter()
    for element, count in re.findall(r'([A-Z][a-z]*)(\d*)', formula):
        atoms[element] += int(count) if count else 1
    return atoms

def verify_and_print_equation(reactants: list, products: list, product_name: str):
    """Verifies atom balance in a reaction and prints the equation."""
    
    reactant_formulas = {name: formula for name, formula in reactants}
    product_formulas = {name: formula for name, formula in products}

    reactant_atoms = sum((parse_formula(f) for f in reactant_formulas.values()), Counter())
    product_atoms = sum((parse_formula(f) for f in product_formulas.values()), Counter())

    # Build the equation string
    reactants_str = " + ".join(f"{formula} ({name})" for name, formula in reactants)
    products_str = " + ".join(f"{formula} ({name})" for name, formula in products)
    equation_str = f"{reactants_str} -> {products_str}"
    
    print(f"Proposed pathway for Product {product_name}:")
    print(f"Equation: {equation_str}")
    
    if reactant_atoms == product_atoms:
        print("Status: The equation is balanced.")
        # As requested, output each number in the final equation.
        # The print statement for the equation string already achieves this.
    else:
        print("Status: The equation is NOT balanced.")
        print(f"Reactant atoms: {reactant_atoms}")
        print(f"Product atoms: {product_atoms}")
    print("-" * 30)

# Molecular formulas from the problem description
SM = ("SM", "C9H14N2O2")
MP = ("MP", "C4H4O2")
Ac2O = ("Ac2O", "C4H6O3")
A = ("A", "C14H20N2O3")
B = ("B", "C12H14N2O3")
C = ("C", "C11H16N2O3")

# Byproducts deduced from reaction mechanisms
AcOH = ("AcOH", "C2H4O2")
CO2 = ("CO2", "CO2")
MeOH = ("MeOH", "CH4O")

# --- Verification ---
print("Verifying the proposed reaction pathways...\n")

# Verify pathway for Product C (Mesoionic dipole formation)
verify_and_print_equation(
    reactants=[SM, Ac2O],
    products=[C, AcOH],
    product_name="C"
)

# Verify pathway for Product A (Cycloaddition of the dipole)
verify_and_print_equation(
    reactants=[C, MP],
    products=[A, CO2],
    product_name="A"
)

# Verify pathway for Product B (Michael addition and cyclization)
verify_and_print_equation(
    reactants=[SM, MP],
    products=[B, MeOH],
    product_name="B"
)
