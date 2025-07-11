# Number of independent surface code patches
num_patches = 2

# Number of holes in each patch
num_holes_per_patch = 2

# The number of logical qubits (k) that can be encoded in a single
# planar surface code patch is equal to the number of holes (h).
# Formula: k = h
logical_qubits_per_patch = num_holes_per_patch

# For multiple independent patches, the total number of logical qubits is the sum
# of the qubits from each patch. Since they are identical, we multiply.
total_logical_qubits = num_patches * logical_qubits_per_patch

print(f"Number of logical qubits per patch with {num_holes_per_patch} holes: {logical_qubits_per_patch}")
print("Total number of logical qubits for two such patches is the sum of the qubits from each.")
print("The final calculation is:")
print(f"{logical_qubits_per_patch} + {logical_qubits_per_patch} = {total_logical_qubits}")