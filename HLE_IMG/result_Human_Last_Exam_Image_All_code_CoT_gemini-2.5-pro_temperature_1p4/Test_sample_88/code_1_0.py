import re
from collections import Counter

def parse_formula(formula):
    """Parses a chemical formula string into a dictionary of atom counts."""
    pairs = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    atom_counts = Counter()
    for element, count in pairs:
        atom_counts[element] += int(count) if count else 1
    return atom_counts

def verify_and_print_reaction(reactants, products, name):
    """
    Verifies the atom balance for a reaction and prints the equation.
    Reactants and products are lists of formula strings.
    """
    reactant_atoms = Counter()
    for formula in reactants:
        reactant_atoms.update(parse_formula(formula))

    product_atoms = Counter()
    for formula in products:
        product_atoms.update(parse_formula(formula))

    # Build the equation string
    reactant_str = " + ".join(reactants)
    product_str = " + ".join(products)
    equation = f"{reactant_str} -> {product_str}"

    print(f"Proposed reaction for {name}:")
    print(equation)
    if reactant_atoms == product_atoms:
        print("Status: The equation is balanced.")
    else:
        print("Status: The equation is NOT balanced.")
        print(f"Reactant Atoms: {dict(reactant_atoms)}")
        print(f"Product Atoms: {dict(product_atoms)}")
    print("-" * 30)

# --- Define Molecular Formulas ---
SM = 'C9H14N2O2'    # Starting Material
MP = 'C4H4O2'       # Methyl Propiolate
AC2O = 'C4H6O3'     # Acetic Anhydride

PRODUCT_A = 'C14H20N2O3'
PRODUCT_B = 'C12H14N2O3'
PRODUCT_C = 'C11H16N2O3'

# Byproducts deduced from reaction mechanisms
ACETIC_ACID = 'C2H4O2'
METHANOL = 'CH4O'
CARBON_DIOXIDE = 'CO2'


# --- Verify Reaction Pathways ---

# 1. Verification for Product C (N-acetylation)
reactants_C = [SM, AC2O]
products_C = [PRODUCT_C, ACETIC_ACID]
verify_and_print_reaction(reactants_C, products_C, "Product C")

# 2. Verification for Product B (Michael Addition + Lactonization)
reactants_B = [SM, MP]
products_B = [PRODUCT_B, METHANOL]
verify_and_print_reaction(reactants_B, products_B, "Product B")

# 3. Verification for Product A ([3+2] Cycloaddition)
# This reaction proceeds from Product C
reactants_A = [PRODUCT_C, MP]
products_A = [PRODUCT_A, CARBON_DIOXIDE]
verify_and_print_reaction(reactants_A, products_A, "Product A")
