import math

# Step 1: Define the constants based on the problem description.
sheet_w = 140  # cm
sheet_h = 110  # cm

square_side = 10  # cm
circle_radius = 20  # cm
circle_bounding_box = circle_radius * 2  # cm

chars_per_square = 4
chars_per_circle = 999

print("--- Solving the Emperor's Engraving Problem ---")
print(f"The goal is to maximize K = (N * {chars_per_square}) + (M * {chars_per_circle})")
print(f"by cutting N squares ({square_side}x{square_side}cm) and M circles (r={circle_radius}cm) from a {sheet_w}x{sheet_h}cm sheet.\n")

# Step 2: Prioritize circles to maximize M, as they hold more characters.
# We will use the circle's bounding box (40x40 cm) for cutting.
print("Step 1: Calculate the maximum number of circles (M).")
circles_along_w = math.floor(sheet_w / circle_bounding_box)
circles_along_h = math.floor(sheet_h / circle_bounding_box)
M = circles_along_w * circles_along_h
print(f"The sheet can fit {circles_along_w} circles along the {sheet_w}cm side and {circles_along_h} circles along the {sheet_h}cm side.")
print(f"Maximum number of circles (M) = {circles_along_w} * {circles_along_h} = {M}\n")


# Step 3: Calculate the area used by circles and the remaining area for squares.
print("Step 2: Calculate the number of squares (N) from the remaining area.")
used_w_by_circles = circles_along_w * circle_bounding_box
used_h_by_circles = circles_along_h * circle_bounding_box
print(f"The {M} circles will be cut from a {used_w_by_circles}x{used_h_by_circles}cm block.")

# The remaining L-shaped area is split into two rectangles for optimal cutting.
# Rectangle 1
rem_rect1_w = sheet_w - used_w_by_circles
rem_rect1_h = sheet_h
# Rectangle 2
rem_rect2_w = used_w_by_circles
rem_rect2_h = sheet_h - used_h_by_circles

print(f"The remaining area is an L-shape, which can be split into two rectangles:")
print(f" - Rectangle 1: {rem_rect1_w}cm x {rem_rect1_h}cm")
print(f" - Rectangle 2: {rem_rect2_w}cm x {rem_rect2_h}cm")

# Calculate how many squares fit in each remaining rectangle.
squares_in_rect1 = math.floor(rem_rect1_w / square_side) * math.floor(rem_rect1_h / square_side)
squares_in_rect2 = math.floor(rem_rect2_w / square_side) * math.floor(rem_rect2_h / square_side)
N = squares_in_rect1 + squares_in_rect2

print(f"Number of squares from Rectangle 1 = {math.floor(rem_rect1_w / square_side)} * {math.floor(rem_rect1_h / square_side)} = {squares_in_rect1}")
print(f"Number of squares from Rectangle 2 = {math.floor(rem_rect2_w / square_side)} * {math.floor(rem_rect2_h / square_side)} = {squares_in_rect2}")
print(f"Maximum number of squares (N) = {squares_in_rect1} + {squares_in_rect2} = {N}\n")

# Step 4: Calculate the total number of characters, K.
print("Step 3: Calculate the total number of characters (K).")
K = (N * chars_per_square) + (M * chars_per_circle)
print(f"K = (N * chars_per_square) + (M * chars_per_circle)")
print(f"K = ({N} * {chars_per_square}) + ({M} * {chars_per_circle})")
print(f"K = {N * chars_per_square} + {M * chars_per_circle} = {K}\n")

# Step 5: Output the final answer.
print("--- Final Answer ---")
print("The optimal production plan is N:M:K")
print(f"{N}:{M}:{K}")
<<<58:6:6226>>>