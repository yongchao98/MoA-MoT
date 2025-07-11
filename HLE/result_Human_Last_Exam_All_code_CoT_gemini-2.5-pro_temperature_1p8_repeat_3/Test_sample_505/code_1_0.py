import sys

# Define the given parameters from the problem description
num_patches = 2
total_holes = 2

# The number of logical qubits, k, that can be encoded in a single
# planar patch of a surface code is equal to the number of holes, h,
# that the patch contains.
# Formula for a single patch: k = h

# For a system with multiple, separate patches, the total number of
# logical qubits (K_total) is the sum of the logical qubits from each patch.
# Let's say we have patch 1 and patch 2.
# K_total = k_patch1 + k_patch2
#
# Using the principle from above, we can substitute k with h for each patch:
# K_total = h_patch1 + h_patch2
#
# The sum of holes in each patch is, by definition, the total number of holes
# in the system.
# Therefore, K_total = total_holes.

total_logical_qubits = total_holes

# Print the reasoning and the final equation.
print(f"The system has {num_patches} patches and a total of {total_holes} holes.")
print("The number of logical qubits on a single surface code patch is equal to its number of holes.")
print("For multiple patches, the total number of logical qubits is the sum of holes across all patches.")
print()
print("Final Equation:")
print(f"Maximum Logical Qubits = Total Number of Holes")
# We explicitly show the numbers in the final calculation.
print(f"Maximum Logical Qubits = {total_logical_qubits}")

# Redirect final answer to the special format, if not in an interactive environment
if not sys.stdout.isatty():
    print(f"<<<{total_logical_qubits}>>>")
