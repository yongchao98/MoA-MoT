# Number of independent surface code patches.
num_patches = 2

# Number of holes within each patch.
num_holes_per_patch = 2

# In a planar surface code, the maximum number of logical qubits that can be
# encoded is equal to the number of holes created in the surface.
# Each hole supports one logical qubit.
logical_qubits_per_patch = num_holes_per_patch

# The total number of logical qubits is the product of the number of patches
# and the number of logical qubits each patch can encode.
total_logical_qubits = num_patches * logical_qubits_per_patch

# Print the final calculation as an equation.
# The equation shows: (number of patches) * (logical qubits per patch) = total logical qubits.
print(f"The final calculation is:")
print(f"{num_patches} patches * {logical_qubits_per_patch} logical qubits/patch = {total_logical_qubits} total logical qubits")