# Define constants from the problem statement
SHEET_W, SHEET_H = 140, 110
CIRCLE_D, SQUARE_S = 40, 10
CIRCLE_VALUE, SQUARE_VALUE = 999, 4

# Step 1: Calculate the maximum number of circles (M)
# A 3x2 grid of 40x40cm squares for the circles fits on the 140x110cm sheet.
num_circles_w = SHEET_W // CIRCLE_D
num_circles_h = SHEET_H // CIRCLE_D
M = num_circles_w * num_circles_h

# Step 2: Calculate the maximum number of squares (N) from the leftover material
# The area used by the block of circles
used_w = num_circles_w * CIRCLE_D
used_h = num_circles_h * CIRCLE_D

# The leftover L-shaped area is decomposed into two rectangles.
# We calculate the number of squares that can be cut from each.
rem_rect1_w = SHEET_W - used_w
rem_rect1_h = SHEET_H
squares_1 = (rem_rect1_w // SQUARE_S) * (rem_rect1_h // SQUARE_S)

rem_rect2_w = used_w
rem_rect2_h = SHEET_H - used_h
squares_2 = (rem_rect2_w // SQUARE_S) * (rem_rect2_h // SQUARE_S)

# The total number of squares is the sum from the two leftover rectangles.
N = squares_1 + squares_2

# Step 3: Calculate the total number of engraved characters (K)
K = (SQUARE_VALUE * N) + (CIRCLE_VALUE * M)

# Print the final results, including the equation with all numbers
print("The optimal production numbers are:")
print(f"N (number of squares) = {N}")
print(f"M (number of circles) = {M}")
print("\nThe equation for the maximum total characters (K) is:")
print(f"{SQUARE_VALUE} * {N} + {CIRCLE_VALUE} * {M} = {K}")
print("\nFinal Answer (N:M:K):")
print(f"{N}:{M}:{K}")