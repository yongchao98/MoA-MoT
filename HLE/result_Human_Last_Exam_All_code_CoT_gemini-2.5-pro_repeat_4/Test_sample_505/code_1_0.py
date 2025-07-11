# Define the parameters of the system
num_patches = 2
total_holes = 2

# The formula for logical qubits 'k' in a single patch with 'h' holes is k = 1 + h.
# The total number of logical qubits is the sum from all patches.
# k_total = k_patch1 + k_patch2 = (1 + h1) + (1 + h2)
# Since h1 + h2 must equal total_holes, the result is constant.
# We can demonstrate this by choosing one possible distribution for the holes.
# Let's put all holes in the first patch.
holes_in_patch1 = 2
holes_in_patch2 = 0

# Calculate the number of logical qubits for each patch
qubits_in_patch1 = 1 + holes_in_patch1
qubits_in_patch2 = 1 + holes_in_patch2

# Calculate the total number of logical qubits
total_logical_qubits = qubits_in_patch1 + qubits_in_patch2

print("To find the maximum number of logical qubits, we calculate the total based on the number of patches and holes.")
print("We have 2 patches and 2 holes in total.")
print("\nLet's use one possible distribution: {} holes in Patch 1 and {} holes in Patch 2.".format(holes_in_patch1, holes_in_patch2))
print("The number of logical qubits for each patch is calculated as '1 + number_of_holes'.")
print("\nThe final equation to calculate the total is:")
print("Total Logical Qubits = (1 + {}) + (1 + {}) = {}".format(holes_in_patch1, holes_in_patch2, total_logical_qubits))
