# Plan: Decompose 2A - 3B as 2(A-B) - B
# We will calculate the cost of this operation sequence, assuming
# inputs A and B are in affine coordinates and the output must be
# in extended coordinates.
# Squaring is assumed to have the same cost as multiplication (1M).

# Step 1: Compute C = A - B.
# This is A_affine + (-B)_affine -> C_projective.
# Cost is 8 multiplications (8M).
cost_add_aff_to_proj = 8

# Step 2: Compute D = 2C.
# This is a doubling in projective coordinates: C_projective -> D_projective.
# Cost is 8 multiplications (8M).
cost_dbl_proj = 8

# Step 3: Convert D from projective to extended coordinates.
# D_projective -> D_extended. A point (X:Y:Z) becomes (XZ:YZ:Z^2:XY).
# This costs 3 multiplications and 1 squaring.
# Total cost = 3M + 1S = 4M.
cost_proj_to_ext = 4

# Step 4: Compute E = D - B.
# This is D_extended + (-B)_affine -> E_extended.
# This is a mixed addition operation.
# Cost is 9 multiplications (9M).
cost_madd_ext_aff = 9

# Calculate the total cost by summing the costs of each step.
total_cost = cost_add_aff_to_proj + cost_dbl_proj + cost_proj_to_ext + cost_madd_ext_aff

# Print the final equation showing the cost of each step.
print(f"The total cost is calculated by summing the cost of each step:")
print(f"{cost_add_aff_to_proj}M (for A-B) + {cost_dbl_proj}M (for 2(A-B)) + {cost_proj_to_ext}M (for conversion) + {cost_madd_ext_aff}M (for final subtraction) = {total_cost}M")
