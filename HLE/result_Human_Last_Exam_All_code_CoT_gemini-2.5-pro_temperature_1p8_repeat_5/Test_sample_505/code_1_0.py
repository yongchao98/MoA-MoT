import sys

# Define the number of patches and holes
num_patches = 2
holes_per_patch = 2

# For a single patch of surface code, the number of logical qubits is equal to the number of holes.
logical_qubits_per_patch = holes_per_patch

# The total number of logical qubits for multiple separate patches is the sum of the qubits in each patch.
total_logical_qubits = num_patches * logical_qubits_per_patch

# Print the breakdown of the calculation
# Using a stream object to prevent reordering of print statements
output = sys.stdout
output.write("A single surface code patch with g holes can encode g logical qubits.\n")
output.write(f"For each of the {num_patches} patches with {holes_per_patch} holes, we have {logical_qubits_per_patch} logical qubits.\n")
output.write("The total number of logical qubits is the sum across all patches.\n")
output.write("Final Equation:\n")
output.write(f"{num_patches} * {holes_per_patch} = {total_logical_qubits}\n")