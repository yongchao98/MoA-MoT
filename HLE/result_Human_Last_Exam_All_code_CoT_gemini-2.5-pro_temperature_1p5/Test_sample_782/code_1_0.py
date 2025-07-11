# Plan to calculate the minimum cost for 2A - 3B on a twisted Edwards curve.
# The cost is measured in the number of multiplications (M) in the underlying field.
# A and B are given in affine coordinates, the result is in extended coordinates.

# --- Operation Costs ---
# Cost to convert an affine point (x, y) to extended coordinates (x, y, 1, x*y).
cost_affine_to_ext = 1 # 1 multiplication for T = x*y

# Cost for mixed-coordinate point addition (extended + affine -> extended).
cost_mADD = 7

# Cost for point doubling in extended coordinates.
cost_DBL = 4

# --- Calculation Strategy: R = 2(A - B) - B ---

# This strategy is chosen because it effectively uses mixed-coordinate operations
# to minimize the total number of multiplications.

print("Calculating the cost for R = 2A - 3B using the formula R = 2(A - B) - B")
print("---")

# Step 1: Calculate T = A - B
# This is computed as A + (-B). -B is computed from B at no cost.
# To add A (affine) and -B (affine), we convert A to extended and use mixed addition.
cost_A_minus_B = cost_affine_to_ext + cost_mADD
print(f"Step 1: Cost to compute T = A - B")
print(f"  - Convert A to extended coordinates: {cost_affine_to_ext}M")
print(f"  - Add -B (affine) using mixed addition: {cost_mADD}M")
print(f"  - Subtotal for Step 1: {cost_A_minus_B}M")
print("---")


# Step 2: Calculate U = 2*T
# T is the result from Step 1 and is in extended coordinates. We need to double it.
cost_2T = cost_DBL
print(f"Step 2: Cost to compute U = 2*T")
print(f"  - Double T (extended) using DBL: {cost_DBL}M")
print(f"  - Subtotal for Step 2: {cost_2T}M")
print("---")


# Step 3: Calculate R = U - B
# This is computed as U + (-B). U is in extended coordinates, -B is in affine.
cost_U_minus_B = cost_mADD
print(f"Step 3: Cost to compute R = U - B")
print(f"  - Add -B (affine) to U (extended) using mixed addition: {cost_mADD}M")
print(f"  - Subtotal for Step 3: {cost_U_minus_B}M")
print("---")


# Final Step: Calculate the total cost
total_cost = cost_A_minus_B + cost_2T + cost_U_minus_B

print("Final Calculation:")
# The prompt requires printing the numbers in the final equation.
# Total Cost = (Cost of A to extended) + (Cost of first mADD) + (Cost of DBL) + (Cost of second mADD)
print(f"Total Cost = {cost_affine_to_ext} + {cost_mADD} + {cost_DBL} + {cost_mADD} = {total_cost} Multiplications")
