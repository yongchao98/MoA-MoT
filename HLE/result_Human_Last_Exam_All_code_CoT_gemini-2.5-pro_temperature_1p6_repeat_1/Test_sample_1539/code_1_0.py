import math

# --- Define constants based on the problem description ---
# Dimensions of the meteorite material
RECT_W = 140
RECT_H = 110

# Properties of the circular plates ("Book of Heaven")
CIRC_RADIUS = 20
CIRC_VALUE = 9999
CIRC_DIAMETER = CIRC_RADIUS * 2

# Properties of the square plates ("Book of Earth")
SQUARE_SIDE = 10
SQUARE_VALUE = 360

# --- Strategy: Maximize the number of high-value circles first ---
# This greedy approach is justified because the value density of circles
# (9999 / (40*40) = 6.25) is much higher than squares (360 / (10*10) = 3.6).

# 1. Calculate the maximum number of circles that can be cut in a grid.
# Check the number of 40cm bounding boxes that fit in the 140x110 rectangle.
num_circles_x = math.floor(RECT_W / CIRC_DIAMETER)
num_circles_y = math.floor(RECT_H / CIRC_DIAMETER)
num_circles = num_circles_x * num_circles_y

# 2. Calculate the area used by the grid of circles.
used_width = num_circles_x * CIRC_DIAMETER
used_height = num_circles_y * CIRC_DIAMETER

# 3. Calculate how many squares can be cut from the leftover "L-shaped" area.
# The L-shape is partitioned into two rectangles for easier calculation.
# Partition 1: A rectangle with the full width of the material and leftover height.
leftover_rect1_w = RECT_W
leftover_rect1_h = RECT_H - used_height
squares_in_rect1 = math.floor(leftover_rect1_w / SQUARE_SIDE) * math.floor(leftover_rect1_h / SQUARE_SIDE)

# Partition 2: A rectangle with the leftover width and the height used by circles.
leftover_rect2_w = RECT_W - used_width
leftover_rect2_h = used_height
squares_in_rect2 = math.floor(leftover_rect2_w / SQUARE_SIDE) * math.floor(leftover_rect2_h / SQUARE_SIDE)

# Total number of squares is the sum from both leftover rectangles.
num_squares = squares_in_rect1 + squares_in_rect2

# 4. Calculate the total number of engraved characters (K).
total_value = (num_circles * CIRC_VALUE) + (num_squares * SQUARE_VALUE)

# --- Output the results ---
print(f"Based on the optimal cutting plan:")
print(f"Number of circular plates (N) = {num_circles}")
print(f"Number of square plates (M) = {num_squares}")
print("\nThe final equation for the maximum number of characters (K) is:")
print(f"{num_circles} * {CIRC_VALUE} + {num_squares} * {SQUARE_VALUE} = {total_value}")

# Final answer in the specified N:M:K format.
final_answer = f"{num_circles}:{num_squares}:{total_value}"
print(f"\n<<<{final_answer}>>>")