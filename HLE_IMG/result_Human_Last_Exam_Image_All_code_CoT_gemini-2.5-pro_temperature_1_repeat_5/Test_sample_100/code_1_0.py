def get_formula_string(atoms_dict):
    """Formats a dictionary of atoms into a molecular formula string."""
    # Standard order for organic compounds: C, H, then alphabetical for the rest
    order = ['C', 'H', 'F', 'N', 'O']
    formula = ""
    for element in order:
        if element in atoms_dict:
            count = atoms_dict[element]
            if count > 0:
                formula += element
                if count > 1:
                    formula += str(count)
    return formula

# Step 1: Determine the molecular formula of the starting material.
# Bicyclic core (C6H5N) + Carbonyl (CO) + CF3 group (CF3) + PMB group (C8H9O)
# Let's count atoms directly from the structure:
# C: 6(core) + 1(C=O) + 1(CF3) + 8(PMB) = 16. Let's recount.
# Core bicyclo[2.2.1] skeleton: C1(CH), N2, C3(C=O), C4(CH), C5(C-CF3), C6(CH), C7(CH2). Total: 6 C, 1 N.
# PMB group: p-CH3O-C6H4-CH2-. Total: 8 C, 9 H, 1 O.
# CF3 group: 1 C, 3 F.
# Total C = 6 (from bicycloazaheptenone) + 1 (from CF3) + 8 (from PMB) = 15
# Total H = 1(C1)+1(C4)+1(C6)+2(C7) + 9(PMB) = 14
# Total F = 3
# Total N = 1
# Total O = 1(C=O) + 1(PMB) = 2
atoms_start = {'C': 15, 'H': 14, 'F': 3, 'N': 1, 'O': 2}
print(f"Starting material formula: {get_formula_string(atoms_start)}")

# Step 2: Reaction 1 -> Intermediate 1 (PMB deprotection with CAN)
# This removes the PMB group (C8H9O) and adds one Hydrogen (H) to the Nitrogen.
atoms_inter_1 = atoms_start.copy()
atoms_inter_1['C'] -= 8
atoms_inter_1['H'] -= 9
atoms_inter_1['O'] -= 1
atoms_inter_1['H'] += 1
print(f"Intermediate 1 formula (after PMB deprotection): {get_formula_string(atoms_inter_1)}")

# Step 3: Reaction 2 -> Intermediate 2 (Hydrogenation with Pd/C, H2)
# This reduces the C=C double bond, adding two Hydrogens (H2).
atoms_inter_2 = atoms_inter_1.copy()
atoms_inter_2['H'] += 2
print(f"Intermediate 2 formula (after hydrogenation): {get_formula_string(atoms_inter_2)}")

# Step 4: Reaction 3 -> Product (Lactam hydrolysis with HCl/H2O)
# This adds one molecule of water (H2O) across the amide bond.
atoms_product = atoms_inter_2.copy()
atoms_product['H'] += 2
atoms_product['O'] += 1
print(f"Final product formula (after hydrolysis): {get_formula_string(atoms_product)}")

print("\n--- Final Product Molecular Formula Breakdown ---")
# The final equation shows the count for each element in the product.
print(f"Number of Carbon (C) atoms: {atoms_product['C']}")
print(f"Number of Hydrogen (H) atoms: {atoms_product['H']}")
print(f"Number of Fluorine (F) atoms: {atoms_product['F']}")
print(f"Number of Nitrogen (N) atoms: {atoms_product['N']}")
print(f"Number of Oxygen (O) atoms: {atoms_product['O']}")