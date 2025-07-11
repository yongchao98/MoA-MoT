# Step 1: Define the number of each structural unit based on visual analysis of the molecule.
num_phenylene_units = 12
num_ethenylene_units = 8
# In a ring, the number of linkers equals the number of units being linked.
num_ethynylene_linkers = num_phenylene_units + num_ethenylene_units

# Step 2: Define the atomic composition of each unit.
# Phenylene unit: C6H4
phenylene_C = 6
phenylene_H = 4
# Ethenylene unit: C2H2
ethenylene_C = 2
ethenylene_H = 2
# Ethynylene linker: C2
ethynylene_C = 2
ethynylene_H = 0

# Step 3: Calculate the total number of Carbon and Hydrogen atoms.
total_C = (num_phenylene_units * phenylene_C) + \
          (num_ethenylene_units * ethenylene_C) + \
          (num_ethynylene_linkers * ethynylene_C)

total_H = (num_phenylene_units * phenylene_H) + \
          (num_ethenylene_units * ethenylene_H) + \
          (num_ethynylene_linkers * ethynylene_H)

# Step 4: Print the analysis and the name of the molecule.
print("Analysis of the molecular structure:")
print(f"Number of phenylene (-C6H4-) units: {num_phenylene_units}")
print(f"Number of ethenylene (-CH=CH-) units: {num_ethenylene_units}")
print(f"Number of ethynylene (-C≡C-) linkers: {num_ethynylene_linkers}")
print("\nCalculating the chemical formula:")
print(f"Carbon atoms = ({num_phenylene_units} * {phenylene_C}) + ({num_ethenylene_units} * {ethenylene_C}) + ({num_ethynylene_linkers} * {ethynylene_C}) = {total_C}")
print(f"Hydrogen atoms = ({num_phenylene_units} * {phenylene_H}) + ({num_ethenylene_units} * {ethenylene_H}) + ({num_ethynylene_linkers} * {ethynylene_H}) = {total_H}")
print(f"\nThe chemical formula is C{total_C}H{total_H}.")

# This molecule is known in the scientific literature.
molecule_name = "DBA-NR4"
print(f"The molecule is known as {molecule_name}, which stands for Dehydrobenzo[12]annulene-Based Nanoring with 4-fold symmetry.")
