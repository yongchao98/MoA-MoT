# Plan: Calculate the cost of 2A - 3B by computing it as 2(A - B) - B.
# The cost is measured in the number of field multiplications (M).
# We assume the cost of squaring (S) is equal to one multiplication (S=M).

# Step 1: Compute C = A - B, where A and B are affine.
# The result C is computed directly into extended coordinates.
# Using formula add-2008-hwcd-3, the cost is 9 multiplications.
cost_affine_add_to_extended = 9

# Step 2: Compute D = 2C, where C is in extended coordinates.
# Using formula dbl-2007-bl, the cost is 4S + 4M.
# Since S=M, the cost is 4M + 4M = 8M.
cost_doubling_extended = 8

# Step 3: Compute E = D - B, where D is extended and B is affine.
# This is a mixed-coordinate addition.
# Using formula madd-2007-bl, the cost is 9 multiplications.
cost_mixed_add = 9

# Total cost is the sum of the costs of the three steps.
total_cost = cost_affine_add_to_extended + cost_doubling_extended + cost_mixed_add

# Output the breakdown of the calculation and the final result.
print("The calculation is performed as 2(A - B) - B.")
print("The cost is broken down into three steps:")
print(f"1. Cost of C = A - B (affine to extended): {cost_affine_add_to_extended}M")
print(f"2. Cost of D = 2C (extended doubling): {cost_doubling_extended}M")
print(f"3. Cost of E = D - B (mixed addition): {cost_mixed_add}M")
print("\nThe final cost equation is:")
print(f"{cost_affine_add_to_extended} + {cost_doubling_extended} + {cost_mixed_add} = {total_cost}")

print(f"\nThe smallest cost is {total_cost} multiplications.")