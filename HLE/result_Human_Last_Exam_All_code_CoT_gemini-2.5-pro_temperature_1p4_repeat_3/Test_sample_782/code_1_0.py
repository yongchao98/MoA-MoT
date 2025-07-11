# Plan: Compute 2A - 3B as 2(A - B) - B, which is 2(A + B') + B' where B' = -B.
# The cost of negation B -> B' is negligible.
# A and B are affine points. The result should be in extended coordinates.
# Squaring cost is the same as multiplication (M). Other operations are free.

# Costs for primitive operations on a general twisted Edwards curve in extended coordinates:
# Cost of adding two affine points with the result in extended coordinates.
cost_affine_plus_affine_to_extended = 12

# Cost of doubling a point in extended coordinates.
cost_doubling_extended = 9

# Cost of mixed addition (extended point + affine point).
# Standard cost is 11M, but we can reuse one intermediate value.
cost_mixed_addition_optimized = 11 - 1

# Total cost is the sum of the steps:
# 1. C = A + B' (affine + affine -> extended)
# 2. D = 2C (doubling)
# 3. R = D + B' (mixed addition with optimization)
total_cost = cost_affine_plus_affine_to_extended + cost_doubling_extended + cost_mixed_addition_optimized

print("To compute 2A - 3B, we can use the equivalent and more efficient form 2(A - B) - B.")
print("Let B' = -B. The operation becomes 2(A + B') + B'.")
print("The calculation proceeds in three steps:")
print(f"1. Cost of adding two affine points (A + B') to get an extended point: {cost_affine_plus_affine_to_extended}M")
print(f"2. Cost of doubling the resulting extended point: {cost_doubling_extended}M")
print(f"3. Cost of mixed addition of the result with B' (with optimization): {cost_mixed_addition_optimized}M")
print("\nThe final equation for the total cost is:")
print(f"{cost_affine_plus_affine_to_extended} + {cost_doubling_extended} + {cost_mixed_addition_optimized} = {total_cost}")

print(f"\nThe smallest cost is {total_cost} multiplications.")