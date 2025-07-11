# Define the costs of the primitive operations in terms of field multiplications (M).
# These costs are based on the Explicit-Formulas Database (EFD) for twisted Edwards curves.
# - Cost of adding a projective point and an affine point (mADD_proj).
cost_mADD_proj = 8
# - Cost of doubling a projective point (DBL_proj).
cost_DBL_proj = 7
# - Cost of converting a projective point to extended coordinates.
cost_conversion_proj_to_ext = 4

# The computation is performed as 2(A - B) - B, which is 2(A + B') + B' where B' = -B.
# Step 1: P = A + B' (affine + affine -> projective)
# This is a mixed addition.
step1_cost = cost_mADD_proj

# Step 2: Q = 2P (projective -> projective)
# This is a doubling.
step2_cost = cost_DBL_proj

# Step 3: R = Q + B' (projective + affine -> projective)
# This is another mixed addition.
step3_cost = cost_mADD_proj

# Step 4: Convert final result R from projective to extended coordinates.
step4_cost = cost_conversion_proj_to_ext

# Calculate the total cost
total_cost = step1_cost + step2_cost + step3_cost + step4_cost

# Print the breakdown of the calculation and the final result.
print("The smallest cost is calculated by summing the costs of the chosen operation sequence:")
print(f"Cost of (A+B'): {step1_cost}M")
print(f"Cost of 2*(A+B'): {step2_cost}M")
print(f"Cost of 2*(A+B') + B': {step3_cost}M")
print(f"Cost of final conversion to extended coordinates: {step4_cost}M")
print("\nFinal calculation:")
print(f"{step1_cost} + {step2_cost} + {step3_cost} + {step4_cost} = {total_cost}")

# The final answer as a single value
# print(f"\nFinal Answer: {total_cost}")