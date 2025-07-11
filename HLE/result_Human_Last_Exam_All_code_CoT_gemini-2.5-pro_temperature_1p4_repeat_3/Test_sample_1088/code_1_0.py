# A script to calculate the total number of primitive Gaussians for a molecule
# using the 6-311G** basis set.

# Step 1: Define the number of primitive Gaussians for each atom type
# based on the 6-311G** definition.
primitives_h = 6  # For Hydrogen
primitives_heavy = 12 # For heavy atoms (e.g., C, N, O)

# Step 2: Define the chemical formula of the molecule.
# We will use Methane (CH4) as an example.
# A Methane molecule has 1 Carbon atom and 4 Hydrogen atoms.
num_carbon = 1
num_hydrogen = 4

# Step 3: Calculate the contribution from each type of atom.
total_primitives_c = num_carbon * primitives_heavy
total_primitives_h = num_hydrogen * primitives_h

# Step 4: Calculate the total number of primitives for the molecule.
total_primitives = total_primitives_c + total_primitives_h

# Step 5: Print the breakdown and the final result.
print("Calculation for the 6-311G** basis set:")
print(f"Primitives per Hydrogen atom: {primitives_h}")
print(f"Primitives per heavy atom (like Carbon): {primitives_heavy}")
print("\nFor a Methane (CH4) molecule:")
print(f"The total number of primitive Gaussians is calculated as:")
print(f"({num_carbon} * {primitives_heavy}) + ({num_hydrogen} * {primitives_h}) = {total_primitives}")