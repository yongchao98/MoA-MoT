import math

# Define constants based on the problem description
MATERIAL_WIDTH = 140
MATERIAL_HEIGHT = 110

CIRCLE_RADIUS = 20
CIRCLE_DIAMETER = CIRCLE_RADIUS * 2
CHARS_PER_CIRCLE = 9999

SQUARE_SIDE = 10
CHARS_PER_SQUARE = 360

# --- Step 1: Maximize the number of circular plates ---
# We use simple grid packing. Let's see how many 40x40cm squares (bounding boxes for circles)
# fit into the 140x110cm material.

# Orientation 1: Material as 140x110
circles_along_width1 = MATERIAL_WIDTH // CIRCLE_DIAMETER
circles_along_height1 = MATERIAL_HEIGHT // CIRCLE_DIAMETER
total_circles1 = circles_along_width1 * circles_along_height1

# Orientation 2: Material as 110x140 (just in case)
circles_along_width2 = MATERIAL_HEIGHT // CIRCLE_DIAMETER
circles_along_height2 = MATERIAL_WIDTH // CIRCLE_DIAMETER
total_circles2 = circles_along_width2 * circles_along_height2

# The maximum number of circles is the same in both orientations
N = max(total_circles1, total_circles2)

# --- Step 2: Calculate remaining material for squared plates ---
# Based on the best fit, we arrange a 3x2 grid of circles.
# (140 // 40 = 3, 110 // 40 = 2)
used_width = (MATERIAL_WIDTH // CIRCLE_DIAMETER) * CIRCLE_DIAMETER
used_height = (MATERIAL_HEIGHT // CIRCLE_DIAMETER) * CIRCLE_DIAMETER

# The remaining area is an L-shape that can be split into two rectangles.
# Rectangle A: Under the main block of circles
rem_rect_A_w = used_width
rem_rect_A_h = MATERIAL_HEIGHT - used_height
squares_from_A = (rem_rect_A_w // SQUARE_SIDE) * (rem_rect_A_h // SQUARE_SIDE)

# Rectangle B: Beside the main block of circles
rem_rect_B_w = MATERIAL_WIDTH - used_width
rem_rect_B_h = MATERIAL_HEIGHT
squares_from_B = (rem_rect_B_w // SQUARE_SIDE) * (rem_rect_B_h // SQUARE_SIDE)

# Total number of squared plates
M = squares_from_A + squares_from_B

# --- Step 3: Calculate the maximal number of characters ---
K = (N * CHARS_PER_CIRCLE) + (M * CHARS_PER_SQUARE)

# --- Step 4: Print the final answer ---
print("Maximizing characters by prioritizing more valuable circular plates:")
print(f"Number of circular plates (N): {MATERIAL_WIDTH // CIRCLE_DIAMETER} x {MATERIAL_HEIGHT // CIRCLE_DIAMETER} = {N}")
print(f"Number of squared plates (M) from remaining material: {squares_from_A} + {squares_from_B} = {M}")
print("\nFinal calculation for total characters (K):")
print(f"({N} circular plates * {CHARS_PER_CIRCLE}) + ({M} squared plates * {CHARS_PER_SQUARE}) = {N * CHARS_PER_CIRCLE} + {M * CHARS_PER_SQUARE} = {K}")

print(f"\nThe optimal configuration is N:M:K")
print(f"{N}:{M}:{K}")

# The final answer in the required format
print(f"\n<<<{N}:{M}:{K}>>>")