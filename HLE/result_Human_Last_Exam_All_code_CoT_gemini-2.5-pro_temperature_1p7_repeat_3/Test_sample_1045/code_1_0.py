# Define the number of atoms in toluene (C7H8)
num_carbon = 7
num_hydrogen = 8

# Define the number of contracted Gaussian functions per atom for the 6-31G basis set
# For Carbon (C), this is 9 (1 for 1s core, 4 for inner valence 2s/2p, 4 for outer valence 2s/2p)
cgf_per_c = 9
# For Hydrogen (H), this is 2 (split-valence 1s)
cgf_per_h = 2

# Calculate the total number of functions contributed by each element
total_from_carbon = num_carbon * cgf_per_c
total_from_hydrogen = num_hydrogen * cgf_per_h

# Calculate the final total number of contracted Gaussian functions
total_cgf = total_from_carbon + total_from_hydrogen

# Print the breakdown of the calculation
print(f"For a 6-31G basis set calculation of toluene (C{num_carbon}H{num_hydrogen}):\n")
print(f"Each Carbon atom contributes {cgf_per_c} contracted Gaussian functions.")
print(f"Each Hydrogen atom contributes {cgf_per_h} contracted Gaussian functions.\n")
print("The total number of functions is calculated as follows:")
print(f"Total = (Number of Carbons * Functions per C) + (Number of Hydrogens * Functions per H)")
print(f"Total = ({num_carbon} * {cgf_per_c}) + ({num_hydrogen} * {cgf_per_h})")
print(f"Total = {total_from_carbon} + {total_from_hydrogen}")
print(f"Total = {total_cgf}\n")
print(f"A 6-31G calculation on toluene uses a total of {total_cgf} contracted Gaussian functions.")