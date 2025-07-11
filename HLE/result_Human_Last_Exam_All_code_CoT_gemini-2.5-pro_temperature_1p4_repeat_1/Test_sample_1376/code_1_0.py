import math

# --- Problem Parameters ---

# Dimensions of the meteorite sheet in cm
sheet_w = 140
sheet_h = 110

# Dimensions of the artifacts in cm
square_w = 10
circle_diameter = 40 # A circle of 20cm radius fits in a 40x40cm square bounding box

# Number of characters per artifact
chars_per_square = 4

# Based on the problem description, we deduce the characters per circle
# 999 symbols per circle / 9 symbols per character = 111
chars_per_circle = 111

# --- Step-by-step Calculation ---

# 1. Maximize the number of circles (M) using simple grid packing
#    We prioritize circles because they yield more characters per area unit.
#    (111 chars / 1600 cm^2) > (4 chars / 100 cm^2)
M_in_width = math.floor(sheet_w / circle_diameter)
M_in_height = math.floor(sheet_h / circle_diameter)
M = M_in_width * M_in_height

# 2. Calculate the area occupied by the circles
circles_area_w = M_in_width * circle_diameter
circles_area_h = M_in_height * circle_diameter

# 3. Calculate the number of squares (N) that can fit in the remaining area.
#    The remaining area is an L-shape, which we can split into two rectangles.
#    Rectangle A: The vertical part of the 'L'
rem_area_A_w = sheet_w - circles_area_w
rem_area_A_h = sheet_h
squares_from_A = math.floor(rem_area_A_w / square_w) * math.floor(rem_area_A_h / square_w)

#    Rectangle B: The horizontal part of the 'L'
rem_area_B_w = circles_area_w # Note: This rectangle is under the circles area
rem_area_B_h = sheet_h - circles_area_h
squares_from_B = math.floor(rem_area_B_w / square_w) * math.floor(rem_area_B_h / square_w)

N = squares_from_A + squares_from_B

# 4. Calculate the total number of characters (K)
chars_from_squares = chars_per_square * N
chars_from_circles = chars_per_circle * M
K = chars_from_squares + chars_from_circles

# --- Final Output ---
print("To maximize the number of engraved characters, the workers should produce:")
print(f"- {N} squares")
print(f"- {M} circles")
print("\nThis combination results in a maximal total of K characters.")
print("The calculation is:")
print(f"K = (Characters from Squares) + (Characters from Circles)")
print(f"K = ({chars_per_square} * {N}) + ({chars_per_circle} * {M})")
print(f"K = {chars_from_squares} + {chars_from_circles}")
print(f"K = {K}\n")

print("The final answer in the format N:M:K is:")
print(f"{N}:{M}:{K}")

# The final answer in the requested format
final_answer = f"{N}:{M}:{K}"
# <<<58:6:898>>>