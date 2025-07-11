# Costs for primitive operations on a twisted Edwards curve, with the result in extended coordinates.
# M = Multiplication, S = Squaring. We assume 1S = 1M.
# These costs are based on the Explicit-Formulas Database (EFD).

# Cost to convert an affine point (x,y) to an extended one (x,y,1,xy).
cost_affine_to_extended = 1  # 1M for the x*y product

# Cost for a point doubling of a point in extended coordinates. (3M + 4S)
cost_doubling = 7

# Cost for a mixed-coordinate subtraction: extended_point - affine_point.
cost_mixed_subtraction = 8

# The minimal cost strategy is to compute 2A - 3B as 2(A-B) - B.
# This breaks down into three main operations after the initial conversion.

# Operation 1: Compute D = A - B.
# This involves converting A to extended coords and then subtracting B (affine).
# Since A is an input, we only pay the cost for its first use.
# Note: For clarity, we'll represent the calculation as a sum of fundamental operations.
cost_op1_conversion = cost_affine_to_extended
cost_op1_sub = cost_mixed_subtraction

# Operation 2: Compute E = 2D = 2(A - B).
# This is a doubling of the result from the previous step, which is already in extended coordinates.
cost_op2_doubling = cost_doubling

# Operation 3: Compute F = E - B = 2(A - B) - B.
# This is another mixed-coordinate subtraction.
cost_op3_sub = cost_mixed_subtraction

# The total cost is the sum of the costs of these sequential operations.
total_cost = cost_op1_conversion + cost_op1_sub + cost_op2_doubling + cost_op3_sub

print("The minimal cost is achieved by computing 2A - 3B as 2(A - B) - B.")
print("The cost breakdown is as follows:")
print(f"1. Convert A to extended (cost: {cost_op1_conversion}M) and subtract affine B (cost: {cost_op1_sub}M) to get (A-B).")
print(f"2. Double the result (cost: {cost_op2_doubling}M) to get 2(A-B).")
print(f"3. Subtract affine B again (cost: {cost_op3_sub}M) to get the final result.")
print("\nThe final cost equation is:")
# The problem asks to print the numbers in the final equation.
print(f"{cost_op1_conversion}M + {cost_op1_sub}M + {cost_op2_doubling}M + {cost_op3_sub}M = {total_cost}M")
