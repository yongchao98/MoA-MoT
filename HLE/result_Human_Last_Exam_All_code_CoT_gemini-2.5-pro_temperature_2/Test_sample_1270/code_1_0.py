# A single benzene ring contains 6 carbon atoms.
carbons_per_benzene = 6

# The structure "mercedesbenzene" is interpreted as 1 central ring and 3 attached rings.
num_central_rings = 1
num_attached_rings = 3

# Calculate the total number of carbons.
# Carbons from the central ring.
carbons_in_central_ring = num_central_rings * carbons_per_benzene
# Carbons from the three attached rings.
carbons_in_attached_rings = num_attached_rings * carbons_per_benzene
# The total number of carbon atoms.
total_carbons = carbons_in_central_ring + carbons_in_attached_rings

# Print the final equation showing each number used in the calculation.
print(f"Based on the structure of 1 central benzene ring and 3 attached benzene rings:")
print(f"Total Carbons = ({num_central_rings} * {carbons_per_benzene}) + ({num_attached_rings} * {carbons_per_benzene}) = {total_carbons}")