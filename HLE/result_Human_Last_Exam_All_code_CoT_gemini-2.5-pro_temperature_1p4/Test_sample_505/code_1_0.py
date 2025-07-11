# Number of separate surface code patches
num_patches = 2

# Number of holes in each patch
holes_per_patch = 2

# For a planar surface code, the number of logical qubits is equal to the number of holes.
qubits_per_patch = holes_per_patch

# Since the patches are independent, the total number of logical qubits is the sum.
total_logical_qubits = num_patches * qubits_per_patch

# Print the final equation
# We will show the addition of qubits from each patch.
patch_contributions = [str(qubits_per_patch) for _ in range(num_patches)]
equation_str = " + ".join(patch_contributions)
print(f"The total number of logical qubits is the sum of qubits from each patch:")
print(f"{equation_str} = {total_logical_qubits}")