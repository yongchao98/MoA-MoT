# Number of carbons in the central benzene ring
carbons_in_central_ring = 6

# Number of attached phenyl groups (representing the Mercedes star)
num_phenyl_groups = 3

# Number of carbons in each phenyl group
carbons_per_phenyl_group = 6

# Calculate the total carbons from the attached groups
total_carbons_in_groups = num_phenyl_groups * carbons_per_phenyl_group

# Calculate the total number of carbons in the molecule
total_carbons = carbons_in_central_ring + total_carbons_in_groups

# Print the equation showing the calculation
# The equation represents: (Carbons in central ring) + (Carbons in 3 phenyl groups) = Total Carbons
print(f"{carbons_in_central_ring} + {total_carbons_in_groups} = {total_carbons}")