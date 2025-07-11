# The product of the Favorskii rearrangement of 2-chlorocyclohexanone is cyclopentanecarboxylic acid.
# The molecular formula is C6H10O2.

# Step 1: Define the counts of each atom.
num_carbon = 6
num_hydrogen = 10
num_oxygen = 2

# Step 2: Calculate the parameters for the Böttcher Molecular Complexity formula.

# N_a: Number of non-hydrogen atoms.
N_a = num_carbon + num_oxygen

# N_b: Number of bonds between non-hydrogen atoms.
# Cyclopentane ring: 5 C-C bonds
# Bond from ring to carboxyl group: 1 C-C bond
# Carboxyl group: 1 C=O bond and 1 C-O bond
# Total bonds between heavy atoms = 5 + 1 + 1 + 1 = 8
N_b = 8

# N_v: Total number of valence electrons.
# Valence electrons: C=4, H=1, O=6
valence_electrons_C = 4
valence_electrons_H = 1
valence_electrons_O = 6
N_v = (num_carbon * valence_electrons_C) + \
      (num_hydrogen * valence_electrons_H) + \
      (num_oxygen * valence_electrons_O)

# Step 3: Calculate the Böttcher Molecular Complexity (BMC).
# Formula: BMC = (N_a * N_b) / (N_v - N_a + 1)
denominator = N_v - N_a + 1
bmc_value = (N_a * N_b) / denominator

# Step 4: Print the results clearly.
print("Calculating the Böttcher Molecular Complexity (BMC) for Cyclopentanecarboxylic Acid (C6H10O2):")
print("-" * 75)
print(f"Number of non-hydrogen atoms (N_a): {N_a}")
print(f"Number of bonds between non-hydrogen atoms (N_b): {N_b}")
print(f"Total number of valence electrons (N_v): {N_v}")
print("\nBMC Formula: (N_a * N_b) / (N_v - N_a + 1)")
print(f"Calculation: ({N_a} * {N_b}) / ({N_v} - {N_a} + 1)")
print(f"             = {N_a * N_b} / {denominator}")
print(f"Final BMC Value: {bmc_value}")
