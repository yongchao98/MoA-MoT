def get_molecular_formula(atoms_dict):
    """Formats a dictionary of atoms into a molecular formula string."""
    formula = ""
    for atom in ['C', 'H', 'O']:
        if atom in atoms_dict and atoms_dict[atom] > 0:
            formula += atom
            if atoms_dict[atom] > 1:
                formula += str(atoms_dict[atom])
    return formula

# Step 1: Define the molecular formula of the starting materials.
# Compound 1 is methyl 1-oxo-2,3-dihydro-1H-indene-2-carboxylate.
# Benzene ring (C6H4) + five-membered ring part (-CO-CH-CH2-) + ester group (-COOCH3)
# C = 6 (aromatic) + 3 (aliphatic ring) + 1 (ester carbonyl) + 1 (methyl) = 11
# H = 4 (aromatic) + 1 (alpha-CH) + 2 (CH2) + 3 (methyl) = 10
# O = 1 (ketone) + 2 (ester) = 3
reactant1_atoms = {'C': 11, 'H': 10, 'O': 3}

# The group added from reactant 2 (benzyl bromide) is a benzyl group.
# Benzyl group = C7H7
benzyl_group = {'C': 7, 'H': 7}

print(f"The molecular formula of Reactant 1 is {get_molecular_formula(reactant1_atoms)}.")
print("The reaction involves two main transformations: alkylation followed by saponification and decarboxylation.")
print("-" * 30)

# Step 2: Calculate the formula of the alkylation intermediate.
# An alpha-hydrogen is replaced by the benzyl group.
# Net change = + C7H7 - H = +C7H6
print("Step A: Alkylation")
print("An alpha-hydrogen is substituted by a benzyl group (C7H7).")
intermediate_atoms = reactant1_atoms.copy()
intermediate_atoms['H'] -= 1  # Remove alpha-H
intermediate_atoms['C'] += benzyl_group['C']  # Add benzyl C
intermediate_atoms['H'] += benzyl_group['H']  # Add benzyl H
print(f"Formula of alkylated intermediate: {get_molecular_formula(intermediate_atoms)}")
print("-" * 30)


# Step 3: Calculate the formula of the final product A.
# The methoxycarbonyl group (-COOMe) is removed and replaced with a hydrogen atom.
# Atoms removed: C2H3O2. Atoms added: H.
print("Step B: Saponification and Decarboxylation")
print("The methoxycarbonyl group (-COOCH3, which is C2H3O2) is replaced by a hydrogen atom (H).")
product_A_atoms = intermediate_atoms.copy()
product_A_atoms['C'] -= 2 # from -COOCH3
product_A_atoms['H'] -= 3 # from -COOCH3
product_A_atoms['O'] -= 2 # from -COOCH3
product_A_atoms['H'] += 1 # add H
print("-" * 30)

# Print the final result
c = product_A_atoms['C']
h = product_A_atoms['H']
o = product_A_atoms['O']
final_formula = get_molecular_formula(product_A_atoms)

print("Final Calculation for Product A:")
print(f"  Carbon atoms = {intermediate_atoms['C']} - 2 = {c}")
print(f"  Hydrogen atoms = {intermediate_atoms['H']} - 3 + 1 = {h}")
print(f"  Oxygen atoms = {intermediate_atoms['O']} - 2 = {o}")
print(f"\nThe molecular formula of compound A is {final_formula}.")
