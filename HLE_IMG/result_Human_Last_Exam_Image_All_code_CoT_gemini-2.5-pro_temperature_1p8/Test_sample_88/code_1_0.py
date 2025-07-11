import re

def parse_formula(formula_str):
    """Parses a molecular formula string into a dictionary of atom counts."""
    atoms = {'C': 0, 'H': 0, 'N': 0, 'O': 0}
    # Find all atom-count pairs (e.g., C10, H16, N2, O2)
    for atom, count in re.findall(r'([A-Z][a-z]?)(\d*)', formula_str):
        atoms[atom] = int(count) if count else 1
    return atoms

def format_formula(atom_dict):
    """Formats a dictionary of atom counts back into a molecular formula string."""
    return ''.join([f"{atom}{count}" for atom, count in atom_dict.items() if count > 0])

def add_formulas(f1, f2):
    """Adds two molecular formulas (represented as dictionaries)."""
    result = f1.copy()
    for atom, count in f2.items():
        result[atom] = result.get(atom, 0) + count
    return result

def subtract_formulas(f1, f2):
    """Subtracts one molecular formula from another."""
    result = f1.copy()
    for atom, count in f2.items():
        result[atom] = result.get(atom, 0) - count
    return result

# --- Given Molecular Formulas ---
SM_formula_str = "C10H16N2O2"  # Proline Dimer
MP_formula_str = "C4H4O2"      # Methyl Propiolate
acetyl_group_str = "C2H3O"     # Acetyl group from Ac2O
H_atom_str = "H1"

product_A_str = "C14H20N2O3"
product_B_str = "C12H14N2O3"
product_C_str = "C11H16N2O3"

# --- Parse all formulas into dictionaries ---
SM = parse_formula(SM_formula_str)
MP = parse_formula(MP_formula_str)
acetyl_group = parse_formula(acetyl_group_str)
H = parse_formula(H_atom_str)

print("--- Verifying formation of Product B (C12H14N2O3) ---")
print(f"Step 1: Acetylation of Starting Material (SM).")
print("An acetyl group (C2H3O) from acetic anhydride replaces an H on a nitrogen atom.")
print(f"Equation: {SM_formula_str} + {acetyl_group_str} - H = N-acetyl-SM")
# N-acetylation: Add acetyl group, remove one H
n_acetyl_sm = subtract_formulas(add_formulas(SM, acetyl_group), H)
print(f"Resulting N-acetyl-SM formula: C10 H16 N2 O2 + C2 H3 O - H = {format_formula(n_acetyl_sm)}")
print("-" * 20)

print(f"Step 2: Oxidation of N-acetyl-SM.")
print("The intermediate loses 4 Hydrogen atoms (2 molecules of H2).")
print(f"Equation: {format_formula(n_acetyl_sm)} - 4*H = Product B")
# Oxidation: Remove 4 H atoms
product_B_calculated = subtract_formulas(n_acetyl_sm, {'C':0, 'H':4, 'N':0, 'O':0})
print(f"Result: {format_formula(n_acetyl_sm)} - H4 = {format_formula(product_B_calculated)}")
print(f"This matches Product B's given formula: {product_B_str}\n")


print("\n--- Verifying formation of Product A (C14H20N2O3) ---")
print(f"Step 1: Michael Addition of SM to Methyl Propiolate (MP).")
print("Equation: SM + MP = Adduct")
adduct_A = add_formulas(SM, MP)
print(f"Result: {SM_formula_str} + {MP_formula_str} = {format_formula(adduct_A)}")
print("-" * 20)

print("Step 2: Reduction of the adduct.")
print("The adduct loses one Oxygen atom.")
print(f"Equation: {format_formula(adduct_A)} - O = Product A")
product_A_calculated = subtract_formulas(adduct_A, {'C':0, 'H':0, 'N':0, 'O':1})
print(f"Result: {format_formula(adduct_A)} - O = {format_formula(product_A_calculated)}")
print(f"This matches Product A's given formula: {product_A_str}\n")


print("\n--- Verifying formation of Product C (C11H16N2O3) ---")
print("This pathway involves the previously calculated N-acetyl-SM intermediate.")
print(f"Step 1: Decarboxylation of N-acetyl-SM ({format_formula(n_acetyl_sm)}).")
print("The intermediate loses a molecule of carbon dioxide (CO2).")
print(f"Equation: {format_formula(n_acetyl_sm)} - CO2 = Intermediate")
intermediate_C = subtract_formulas(n_acetyl_sm, parse_formula("CO2"))
print(f"Result: {format_formula(n_acetyl_sm)} - CO2 = {format_formula(intermediate_C)}")
print("-" * 20)

print("Step 2: Oxidation of the intermediate.")
print("The intermediate is oxidized, losing 2 Hydrogen atoms and gaining 2 Oxygen atoms.")
print(f"Equation: {format_formula(intermediate_C)} - 2*H + 2*O = Product C")
product_C_calculated = add_formulas(subtract_formulas(intermediate_C, {'C':0, 'H':2, 'N':0, 'O':0}), {'C':0, 'H':0, 'N':0, 'O':2})
print(f"Result: {format_formula(intermediate_C)} - H2 + O2 = {format_formula(product_C_calculated)}")
print(f"This matches Product C's given formula: {product_C_str}\n")
