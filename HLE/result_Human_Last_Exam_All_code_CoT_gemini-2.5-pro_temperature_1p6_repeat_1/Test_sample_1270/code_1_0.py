# A benzene ring has 6 carbon atoms.
carbons_per_ring = 6

# "mercedesbenzene" is interpreted as a central benzene ring.
carbons_in_central_ring = carbons_per_ring

# It has three other benzene rings attached, like the Mercedes logo.
num_attached_rings = 3

# Calculate the total number of carbons.
total_carbons = carbons_in_central_ring + (num_attached_rings * carbons_per_ring)

# Print the final equation showing how the total is calculated.
print(f"A central benzene ring plus three attached rings gives a total of:")
print(f"{carbons_in_central_ring} + {num_attached_rings} * {carbons_per_ring} = {total_carbons} carbons")