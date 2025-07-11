# The task is to find the minimum cost to compute 2A - 3B on a twisted Edwards curve.
#
# Assumptions:
# - Input points A and B are in affine coordinates.
# - The final result is required in extended coordinates.
# - The cost of a squaring (S) is the same as a multiplication (M).
# - The curve parameter a = -1.
# - Costs are taken from standard benchmarks for elliptic curve cryptography (e.g., EFD).
#
# Strategy:
# The computation is performed using the identity 2A - 3B = 2(A-B) - B, which is
# generally more efficient than computing 2A and 3B separately.
#
# Cost Breakdown:
# 1. Cost of D = A - B: To add two affine points and get an extended result, we
#    convert each to extended coordinates (1M each for the T coordinate) and then
#    add them (8M).
# 2. Cost of 2D: Doubling a point in extended coordinates costs 3M + 4S = 7M.
# 3. Cost of 2D - B: This is a mixed-mode addition (extended + affine), which costs 9M.

# Cost of adding two affine points (A and -B) to get an extended result D.
# Cost is 1M (for A's T-coord) + 1M (for -B's T-coord) + 8M (for extended addition).
cost_A_minus_B = 1 + 1 + 8

# Cost of doubling the point D, which is in extended coordinates.
# This operation costs 3M + 4S. Assuming S=M, cost is 7M.
cost_doubling = 7

# Cost of the final mixed addition (2D - B), where 2D is extended and B is affine.
# This operation costs 9M.
cost_mixed_add = 9

# Calculate the total cost by summing the costs of each step.
total_cost = cost_A_minus_B + cost_doubling + cost_mixed_add

print("To compute 2A - 3B with minimum cost, we use the identity 2(A-B) - B.")
print("The cost breakdown is as follows:")
print(f"1. Cost of (A - B) from affine to extended: {cost_A_minus_B}M")
print(f"2. Cost of doubling the result: {cost_doubling}M")
print(f"3. Cost of final mixed-addition subtraction: {cost_mixed_add}M")
print(f"Final equation for the total cost: {cost_A_minus_B} + {cost_doubling} + {cost_mixed_add} = {total_cost}")
print(f"The smallest cost is {total_cost} multiplications.")
