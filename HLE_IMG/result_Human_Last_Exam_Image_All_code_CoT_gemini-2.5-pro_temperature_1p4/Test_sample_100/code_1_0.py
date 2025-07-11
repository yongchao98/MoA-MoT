def format_formula(formula_dict):
    """Formats a dictionary of atoms into a standard molecular formula string."""
    # Standard order: C, H, then alphabetical for the rest
    order = ['C', 'H', 'F', 'N', 'O']
    parts = []
    for atom in order:
        count = formula_dict.get(atom, 0)
        if count > 0:
            parts.append(f"{atom}{count if count > 1 else ''}")
    return "".join(parts)

def apply_reaction(current_formula, atoms_to_remove, atoms_to_add):
    """Applies the changes of a reaction to a molecular formula."""
    next_formula = current_formula.copy()
    for atom, count in atoms_to_remove.items():
        if atom in next_formula:
            next_formula[atom] -= count
    for atom, count in atoms_to_add.items():
        next_formula[atom] = next_formula.get(atom, 0) + count
    return next_formula

# --- Step 0: Define the starting material ---
# Core (2-azabicyclo[2.2.1]hept-5-en-3-one): C6H7NO
# Substituents: PMB (C8H9O) on N, CF3 (CF3) on C=C.
# Substitution removes one H from N and one H from C.
# Total C = 6(core) + 8(PMB) + 1(CF3) = 15
# Total H = 7(core) + 9(PMB) - 1(for PMB sub) - 1(for CF3 sub) = 14
# Total F = 3
# Total N = 1
# Total O = 1(core) + 1(PMB) = 2
starting_material = {'C': 15, 'H': 14, 'F': 3, 'N': 1, 'O': 2}
print(f"Step 0: The determined starting material formula is {format_formula(starting_material)}.")

# --- Step 1: Deprotection of PMB group ---
pmb_group = {'C': 8, 'H': 9, 'O': 1}
h_atom = {'H': 1}
intermediate_1 = apply_reaction(starting_material, pmb_group, h_atom)
print(f"Step 1 (CAN): Remove PMB group (-C8H9O), Add H (+H). The formula becomes {format_formula(intermediate_1)}.")

# --- Step 2: Hydrogenation ---
h2_molecule = {'H': 2}
intermediate_2 = apply_reaction(intermediate_1, {}, h2_molecule)
print(f"Step 2 (Pd/C, H2): Reduce C=C bond by adding H2 (+H2). The formula becomes {format_formula(intermediate_2)}.")

# --- Step 3: Lactam Hydrolysis ---
h2o_molecule = {'H': 2, 'O': 1}
product = apply_reaction(intermediate_2, {}, h2o_molecule)
print(f"Step 3 (HCl): Hydrolyze the lactam by adding H2O (+H2O). The formula becomes {format_formula(product)}.")

# --- Final Result ---
c = product.get('C', 0)
h = product.get('H', 0)
f = product.get('F', 0)
n = product.get('N', 0)
o = product.get('O', 0)

print("\n--- Final Product Composition ---")
print(f"Number of Carbon (C) atoms: {c}")
print(f"Number of Hydrogen (H) atoms: {h}")
print(f"Number of Fluorine (F) atoms: {f}")
print(f"Number of Nitrogen (N) atoms: {n}")
print(f"Number of Oxygen (O) atoms: {o}")

final_formula = format_formula(product)
print(f"\nThe final molecular formula is: {final_formula}")