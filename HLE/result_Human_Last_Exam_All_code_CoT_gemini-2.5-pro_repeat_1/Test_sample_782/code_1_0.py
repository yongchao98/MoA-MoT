# Define the costs of primitive operations in terms of field multiplications (M).
# Squaring is assumed to have the same cost as multiplication.
# Multiplication by curve constants is assumed to be free.

# Cost to convert a point from affine (x,y) to extended (X:Y:Z:T) coordinates.
# This requires calculating T = x*y.
COST_AFF_TO_EXT = 1

# Cost for mixed-coordinate addition: adding a point in extended coordinates
# to a point in affine coordinates. The result is in extended coordinates.
COST_MADD = 7

# Cost for doubling a point that is already in extended coordinates.
COST_DBL_EXT = 8

print("Calculating the minimum cost for 2A - 3B using the strategy 2(A-B) - B.\n")

# --- Step 1: Calculate P = A - B ---
# This is equivalent to A + (-B). Both A and -B are in affine coordinates.
# To use the efficient mixed addition, we first convert A to extended coordinates.
# Then we perform a mixed addition of A_ext and -B.
# Negating B in affine coordinates, (-x, y), has no multiplication cost.
cost_step1 = COST_AFF_TO_EXT + COST_MADD
print(f"Step 1: Compute P = A - B.")
print(f"  - Convert A to extended coordinates: {COST_AFF_TO_EXT}M")
print(f"  - Perform mixed addition of A_ext and -B: {COST_MADD}M")
print(f"  - Cost of Step 1: {COST_AFF_TO_EXT} + {COST_MADD} = {cost_step1}M\n")

# --- Step 2: Calculate Q = 2P ---
# The result from Step 1, P, is already in extended coordinates.
# We perform a standard doubling operation on an extended point.
cost_step2 = COST_DBL_EXT
print(f"Step 2: Compute Q = 2P.")
print(f"  - Double the extended point P: {COST_DBL_EXT}M")
print(f"  - Cost of Step 2: {cost_step2}M\n")

# --- Step 3: Calculate R = Q - B ---
# This is Q + (-B). Q is in extended coordinates and -B is in affine coordinates.
# This is a mixed addition.
cost_step3 = COST_MADD
print(f"Step 3: Compute final result R = Q - B.")
print(f"  - Perform mixed addition of Q and -B: {COST_MADD}M")
print(f"  - Cost of Step 3: {cost_step3}M\n")

# --- Final Calculation ---
total_cost = cost_step1 + cost_step2 + cost_step3
print("--------------------------------------------------")
print("Total Cost Calculation:")
print(f"Total Cost = (Cost of Step 1) + (Cost of Step 2) + (Cost of Step 3)")
print(f"Total Cost = {cost_step1} + {cost_step2} + {cost_step3} = {total_cost} multiplications.")
print("--------------------------------------------------")