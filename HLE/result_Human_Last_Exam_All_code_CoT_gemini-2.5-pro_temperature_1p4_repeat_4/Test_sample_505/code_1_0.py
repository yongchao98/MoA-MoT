# Number of separate surface code patches
num_patches = 2

# Number of holes in each patch
holes_per_patch = 2

# For a planar surface code, the number of logical qubits is equal to the number of holes.
qubits_per_patch = holes_per_patch

# The total number of logical qubits is the sum of qubits from each independent patch.
total_logical_qubits = qubits_per_patch * num_patches

# Print the breakdown and the final equation.
print(f"A single surface code patch with {holes_per_patch} holes can encode {qubits_per_patch} logical qubits.")
print(f"With {num_patches} such patches, the total number of logical qubits is the sum of the qubits from each patch.")
print(f"Final Equation: {qubits_per_patch} + {qubits_per_patch} = {total_logical_qubits}")