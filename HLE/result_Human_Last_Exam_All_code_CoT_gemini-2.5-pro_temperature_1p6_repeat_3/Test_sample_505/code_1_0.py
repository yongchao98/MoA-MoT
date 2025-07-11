# Number of separate surface code patches
num_patches = 2

# Number of holes in each patch
num_holes_per_patch = 2

# For a single planar surface code patch, the number of logical qubits (k)
# is given by the formula k = N - 1, where N is the total number of boundaries.
# The total number of boundaries is the sum of the outer boundary (1) and the inner boundaries (holes).
num_boundaries_per_patch = 1 + num_holes_per_patch

# Calculate the number of logical qubits for a single patch using the formula k = N - 1
qubits_per_patch = num_boundaries_per_patch - 1

# The total number of logical qubits for the entire system is the sum of the qubits
# from each independent patch. Since the patches are identical, we multiply.
total_logical_qubits = num_patches * qubits_per_patch

# Print the calculation steps and the final result
print(f"Number of boundaries per patch (N) = 1 (outer) + {num_holes_per_patch} (holes) = {num_boundaries_per_patch}")
print(f"Logical qubits per patch = N - 1 = {num_boundaries_per_patch} - 1 = {qubits_per_patch}")
print(f"Total logical qubits = Number of patches * Logical qubits per patch")
print(f"Final Equation: {num_patches} * ((1 + {num_holes_per_patch}) - 1) = {total_logical_qubits}")
print(f"\nTherefore, at most {total_logical_qubits} logical qubits can be encoded.")
