def calculate_dou(c, h, n, x=0):
    """Calculates the Degrees of Unsaturation for a molecule."""
    # Formula: DoU = C + 1 - (H/2) - (X/2) + (N/2)
    return c + 1 - (h / 2) - (x / 2) + (n / 2)

def get_atom_counts(formula):
    """Parses a molecular formula string (e.g., 'C6H11NO3') into a dictionary of atom counts."""
    import re
    atoms = {'C': 0, 'H': 0, 'N': 0, 'O': 0}
    for element, count in re.findall(r'([A-Z])(\d*)', formula):
        atoms[element] = int(count) if count else 1
    return atoms

def print_reaction_balance(reaction_name, reactants, products):
    """Checks and prints the atom balance for a chemical reaction."""
    print(f"--- Verifying Stoichiometry for {reaction_name} ---")

    reactant_atoms = {'C': 0, 'H': 0, 'N': 0, 'O': 0}
    for formula in reactants:
        counts = get_atom_counts(formula)
        for atom, count in counts.items():
            reactant_atoms[atom] += count

    product_atoms = {'C': 0, 'H': 0, 'N': 0, 'O': 0}
    for formula in products:
        counts = get_atom_counts(formula)
        for atom, count in counts.items():
            product_atoms[atom] += count

    print("Total atoms in reactants: ", reactant_atoms)
    print("Total atoms in products:  ", product_atoms)

    if reactant_atoms == product_atoms:
        print("Result: The reaction is balanced.")
    else:
        print("Result: The reaction is NOT balanced.")
    print("-" * 40)


# --- Data from the problem ---
sm1_formula = "C6H11NO3" # N-acetyl-N-methyl-alanine
sm3_formula = "C7H11NO3" # N-acetyl-proline
alkyne_formula = "C4H4O2"  # Methyl propiolate

product_a_formula = "C9H13NO2"
product_b_formula = "C10H13NO2"

byproducts_formula = ["CO2", "H2O"]

# --- Calculations ---

# 1. Calculate Degrees of Unsaturation (DoU) for Products
pa_atoms = get_atom_counts(product_a_formula)
pb_atoms = get_atom_counts(product_b_formula)

dou_a = calculate_dou(pa_atoms['C'], pa_atoms['H'], pa_atoms['N'])
dou_b = calculate_dou(pb_atoms['C'], pb_atoms['H'], pb_atoms['N'])

print("--- Calculating Degrees of Unsaturation (DoU) ---")
print(f"Product A ({product_a_formula}) has a DoU of: {dou_a}")
print(f"Product B ({product_b_formula}) has a DoU of: {dou_b}")
print("-" * 40)


# 2. Verify Atom Balance for the Reactions
# Reaction 1: SM1 + Alkyne -> Product A + CO2 + H2O
print_reaction_balance(
    "Reaction 1",
    reactants=[sm1_formula, alkyne_formula],
    products=[product_a_formula] + byproducts_formula
)

# Reaction 3: SM3 + Alkyne -> Product B + CO2 + H2O
print_reaction_balance(
    "Reaction 3",
    reactants=[sm3_formula, alkyne_formula],
    products=[product_b_formula] + byproducts_formula
)