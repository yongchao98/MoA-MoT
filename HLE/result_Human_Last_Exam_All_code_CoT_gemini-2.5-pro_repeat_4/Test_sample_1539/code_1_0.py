import math

# Define the constants from the problem description
RECT_W = 140
RECT_H = 110
CIRCLE_RADIUS = 20
SQUARE_SIDE = 10

CHARS_PER_CIRCLE = 9999
CHARS_PER_SQUARE = 360

# The bounding box for a circle is a square with a side equal to the circle's diameter
CIRCLE_DIAMETER = CIRCLE_RADIUS * 2

# --- Step 1: Maximize the number of circular plates ---
# We use integer division to find how many 40x40cm squares fit into the 140x110cm rectangle.
# This approach works for both orientations (140x110 and 110x140), yielding the same number of circles.
num_circles_w = RECT_W // CIRCLE_DIAMETER
num_circles_h = RECT_H // CIRCLE_DIAMETER
N = num_circles_w * num_circles_h

# --- Step 2: Calculate remaining area and maximize squared plates ---
# The area used by the circular plates' bounding boxes forms a rectangle.
used_w = num_circles_w * CIRCLE_DIAMETER
used_h = num_circles_h * CIRCLE_DIAMETER

# The remaining material forms an L-shape, which can be split into two rectangles.
# Rectangle 1:
rem_rect1_w = RECT_W - used_w
rem_rect1_h = RECT_H
num_squares1 = (rem_rect1_w // SQUARE_SIDE) * (rem_rect1_h // SQUARE_SIDE)

# Rectangle 2:
rem_rect2_w = used_w
rem_rect2_h = RECT_H - used_h
num_squares2 = (rem_rect2_w // SQUARE_SIDE) * (rem_rect2_h // SQUARE_SIDE)

# Total number of squared plates
M = num_squares1 + num_squares2

# --- Step 3: Calculate the total number of characters ---
K = (N * CHARS_PER_CIRCLE) + (M * CHARS_PER_SQUARE)

# --- Step 4: Print the results as requested ---
print(f"Based on the greedy strategy of maximizing the most valuable plates first:")
print(f"Number of circular plates (N): {N}")
print(f"Number of squared plates (M): {M}")
print("\nThe maximal number of characters (K) is calculated as follows:")
print(f"K = (Number of circles * Chars per circle) + (Number of squares * Chars per square)")
print(f"K = ({N} * {CHARS_PER_CIRCLE}) + ({M} * {CHARS_PER_SQUARE})")
print(f"K = {N * CHARS_PER_CIRCLE} + {M * CHARS_PER_SQUARE}")
print(f"K = {K}")
print("\n--------------------------------------------------")
print(f"The final answer in N:M:K format is: {N}:{M}:{K}")
