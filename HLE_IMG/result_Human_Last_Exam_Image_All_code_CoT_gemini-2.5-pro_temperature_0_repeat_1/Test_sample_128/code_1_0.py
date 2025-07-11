# Define atomic weights for calculation
atomic_weights = {
    'C': 12.011,
    'H': 1.008,
    'N': 14.007,
    'O': 15.999,
}

# Define molecular formulas for the main reactants and the product
reactants_dict = {
    "2-aminopyridine": {'C': 5, 'H': 6, 'N': 2},
    "o-phthalaldehyde": {'C': 8, 'H': 6, 'O': 2},
}
product_A_formula = {'C': 14, 'H': 11, 'N': 3}

# Function to calculate molecular weight
def calculate_mw(formula):
    """Calculates the molecular weight from a dictionary of atom counts."""
    mw = 0
    for atom, count in formula.items():
        mw += atomic_weights.get(atom, 0) * count
    return mw

# Identify Compound A
product_A_name = "2-(pyridin-2-yl)isoindoline-1-carbonitrile"

print(f"The product of the reaction, Compound A, is identified as {product_A_name}.")

# Calculate and print formula and molecular weight of Compound A
formula_str_A = f"C{product_A_formula['C']}H{product_A_formula['H']}N{product_A_formula['N']}"
mw_A = calculate_mw(product_A_formula)
print(f"Molecular Formula: {formula_str_A}")
print(f"Molecular Weight: {mw_A:.3f} g/mol")

# Print the reaction scheme and the atom counts ("numbers") for each species
print("\nThe reaction scheme is:")
r1_formula = f"C{reactants_dict['2-aminopyridine']['C']}H{reactants_dict['2-aminopyridine']['H']}N{reactants_dict['2-aminopyridine']['N']}"
r2_formula = f"C{reactants_dict['o-phthalaldehyde']['C']}H{reactants_dict['o-phthalaldehyde']['H']}O{reactants_dict['o-phthalaldehyde']['O']}"
print(f"{r1_formula} (2-aminopyridine) + {r2_formula} (o-phthalaldehyde) + TMSCN  ->  {formula_str_A} (Compound A)")

print("\nThe numbers (atom counts) for each element in the main chemical species are:")

print("\nReactant: 2-aminopyridine")
for atom, count in reactants_dict["2-aminopyridine"].items():
    print(f"  Number of {atom} atoms: {count}")

print("\nReactant: o-phthalaldehyde")
for atom, count in reactants_dict["o-phthalaldehyde"].items():
    print(f"  Number of {atom} atoms: {count}")

print("\nProduct: Compound A")
for atom, count in product_A_formula.items():
    print(f"  Number of {atom} atoms: {count}")