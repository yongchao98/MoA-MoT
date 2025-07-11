import math

# --- 1. Define constants ---
# Dimensions of the raw material
material_width = 140
material_height = 110

# Dimensions and character counts for the artifacts
square_size = 10
square_chars = 4

circle_radius = 20
circle_diameter = circle_radius * 2
circle_chars = 999

# --- 2. Strategy: Maximize circles first due to higher value density ---
# The bounding box for a circle is a square with side length equal to the diameter.
circle_bbox_size = circle_diameter

# --- 3. Calculate the maximum number of circles (M) ---
# Calculate how many 40x40 bounding boxes fit into the 140x110 material.
# We check both orientations for the grid of circles.

# Orientation 1:
circles_fit_x1 = math.floor(material_width / circle_bbox_size)
circles_fit_y1 = math.floor(material_height / circle_bbox_size)
total_circles_1 = circles_fit_x1 * circles_fit_y1

# Orientation 2 (swapping material dimensions):
circles_fit_x2 = math.floor(material_height / circle_bbox_size)
circles_fit_y2 = math.floor(material_width / circle_bbox_size)
total_circles_2 = circles_fit_x2 * circles_fit_y2

# The maximum number of circles is M
M = max(total_circles_1, total_circles_2)

# Based on the max, determine the grid dimensions
if total_circles_1 > total_circles_2:
    circles_grid_x = circles_fit_x1
    circles_grid_y = circles_fit_y1
else:
    # Handles both total_circles_2 >= total_circles_1
    # For this problem, 3x2 and 2x3 both yield 6. We'll proceed with 3x2.
    circles_grid_x = circles_fit_x1
    circles_grid_y = circles_fit_y1


# --- 4. Calculate remaining area and the number of squares (N) ---
# Space occupied by the grid of circles
occupied_width = circles_grid_x * circle_bbox_size
occupied_height = circles_grid_y * circle_bbox_size

# The remaining L-shaped area can be split into two rectangles.
# Rectangle A: The strip along the width
rem_rect_A_w = material_width - occupied_width
rem_rect_A_h = material_height
squares_in_A = math.floor(rem_rect_A_w / square_size) * math.floor(rem_rect_A_h / square_size)

# Rectangle B: The remaining strip along the height
rem_rect_B_w = occupied_width
rem_rect_B_h = material_height - occupied_height
squares_in_B = math.floor(rem_rect_B_w / square_size) * math.floor(rem_rect_B_h / square_size)

# Total number of squares is N
N = squares_in_A + squares_in_B

# --- 5. Calculate the maximal number of characters (K) ---
K = (N * square_chars) + (M * circle_chars)

# --- 6. Print the final answer ---
print("Optimal Production Plan:")
print(f"Number of squares (N): {N}")
print(f"Number of circles (M): {M}")
print("\nMaximizing the total number of engraved characters (K):")
print(f"K = (N * chars_per_square) + (M * chars_per_circle)")
print(f"K = ({N} * {square_chars}) + ({M} * {circle_chars})")
print(f"K = {N * square_chars} + {M * circle_chars}")
print(f"K = {K}")

print(f"\nFinal Answer in N:M:K format is {N}:{M}:{K}")

<<<58:6:6226>>>