# Number of distinct surface code patches
num_patches = 2

# Number of holes in each patch
num_holes_per_patch = 2

# The formula for logical qubits (k) in a single planar patch is k = h - 1,
# where h is the number of holes (or boundaries).
qubits_per_patch = num_holes_per_patch - 1

# Total number of logical qubits is the sum of qubits from each patch.
total_logical_qubits = num_patches * qubits_per_patch

print(f"For one patch with {num_holes_per_patch} holes, the number of logical qubits is {num_holes_per_patch} - 1 = {qubits_per_patch}.")
print(f"For {num_patches} such patches, the total number of logical qubits is calculated as:")
print(f"{num_patches} * ({num_holes_per_patch} - 1) = {total_logical_qubits}")