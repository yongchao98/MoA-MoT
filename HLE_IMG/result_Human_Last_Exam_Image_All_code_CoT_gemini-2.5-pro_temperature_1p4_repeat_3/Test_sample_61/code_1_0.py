# Step 1: Define the atom counts for the building blocks of the final product A.
# The final product is 2-benzyl-1-indanone.
# We can think of it as an indan-1-one core where one hydrogen is substituted by a benzyl group.

# Molecular formula of indan-1-one is C9H8O.
indan_one_core = {'C': 9, 'H': 8, 'O': 1}

# Molecular formula of a benzyl group is C7H7.
benzyl_group = {'C': 7, 'H': 7, 'O': 0}

# The reaction replaces one hydrogen atom from the indan-1-one core with the benzyl group.
hydrogen_atom = {'C': 0, 'H': 1, 'O': 0}

# Step 2: Calculate the number of atoms for each element in product A.

# Calculate Carbon atoms
c_final = indan_one_core['C'] + benzyl_group['C']
print("Calculation for Carbon (C):")
print(f"Number of C atoms in indan-1-one core = {indan_one_core['C']}")
print(f"Number of C atoms in benzyl group = {benzyl_group['C']}")
print(f"Total C atoms = {indan_one_core['C']} + {benzyl_group['C']} = {c_final}")
print("-" * 20)

# Calculate Hydrogen atoms
h_final = indan_one_core['H'] - hydrogen_atom['H'] + benzyl_group['H']
print("Calculation for Hydrogen (H):")
print(f"Number of H atoms in indan-1-one core = {indan_one_core['H']}")
print(f"Number of H atoms removed = {hydrogen_atom['H']}")
print(f"Number of H atoms in added benzyl group = {benzyl_group['H']}")
print(f"Total H atoms = {indan_one_core['H']} - {hydrogen_atom['H']} + {benzyl_group['H']} = {h_final}")
print("-" * 20)

# Calculate Oxygen atoms
o_final = indan_one_core['O'] + benzyl_group['O']
print("Calculation for Oxygen (O):")
print(f"Number of O atoms in indan-1-one core = {indan_one_core['O']}")
print(f"Number of O atoms in added benzyl group = {benzyl_group['O']}")
print(f"Total O atoms = {indan_one_core['O']} + {benzyl_group['O']} = {o_final}")
print("-" * 20)

# Step 3: Print the final molecular formula.
# We will not display '1' for a single atom.
c_str = f"C{c_final}"
h_str = f"H{h_final}"
o_str = "O" if o_final == 1 else f"O{o_final}"

# Handle the case of zero oxygen atoms
if o_final == 0:
    molecular_formula = f"{c_str}{h_str}"
else:
    molecular_formula = f"{c_str}{h_str}{o_str}"

print(f"The molecular formula of compound A is: {molecular_formula}")