# Number of independent surface code patches
num_patches = 2

# Number of holes in each patch
num_holes_per_patch = 2

# In a planar surface code, the number of logical qubits (k) that can be
# encoded in a single patch is equal to the number of holes (h).
# Formula for one patch: k = h
logical_qubits_per_patch = num_holes_per_patch

# Since the patches are independent, the total number of logical qubits is the
# sum of the qubits from each patch.
total_logical_qubits = logical_qubits_per_patch * num_patches

# Print the final equation showing all the numbers
print(f"Total logical qubits = {logical_qubits_per_patch} logical qubits/patch * {num_patches} patches = {total_logical_qubits}")