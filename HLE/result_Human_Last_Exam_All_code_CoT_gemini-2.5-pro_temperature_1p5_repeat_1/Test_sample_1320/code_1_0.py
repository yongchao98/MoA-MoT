# This script calculates the predicted helix type for an alpha/epsilon-hybrid peptide
# based on an observed trend from similar foldamers.

# Step 1: Define the known helix for a closely related peptide.
# An alpha/delta-peptide forms a helix with two ring sizes, 15 and 17.
alpha_delta_helix_ring1 = 15
alpha_delta_helix_ring2 = 17

# Step 2: Define the change observed when one carbon atom is added to the
# non-alpha-amino acid backbone (e.g., moving from alpha/gamma to alpha/delta).
# The ring sizes increase by 3.
increment = 3

# Step 3: Calculate the predicted ring sizes for the alpha/epsilon-peptide,
# which has one more backbone atom than the alpha/delta-peptide.
predicted_ring1 = alpha_delta_helix_ring1 + increment
predicted_ring2 = alpha_delta_helix_ring2 + increment

# Step 4: Print the final equation showing each number involved in the calculation.
print(f"The calculation for the first ring size is: {alpha_delta_helix_ring1} + {increment} = {predicted_ring1}")
print(f"The calculation for the second ring size is: {alpha_delta_helix_ring2} + {increment} = {predicted_ring2}")
print(f"\nThe predicted helix type is {predicted_ring1}/{predicted_ring2}.")