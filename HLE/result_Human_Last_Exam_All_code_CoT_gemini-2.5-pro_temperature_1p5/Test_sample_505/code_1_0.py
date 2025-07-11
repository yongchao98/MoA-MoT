# Number of separate surface code patches
num_patches = 2

# Number of holes created in each patch
holes_per_patch = 2

# For a planar surface code, the number of logical qubits is equal to the number of holes.
# Therefore, the number of qubits for a single patch is:
qubits_per_patch = holes_per_patch

# The total number of logical qubits is the sum of the qubits from all separate patches.
total_logical_qubits = num_patches * qubits_per_patch

print(f"Each of the {num_patches} patches has {holes_per_patch} holes, allowing each to encode {qubits_per_patch} logical qubits.")
print("The total number of logical qubits is the sum from both patches.")
print(f"Total logical qubits = {qubits_per_patch} + {qubits_per_patch} = {total_logical_qubits}")