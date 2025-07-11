import math

# --- Step 1: Define problem constants ---

# Material dimensions (cm)
material_width_cm = 140
material_height_cm = 110

# "Book of Heaven" - Circle plate properties
circle_radius_cm = 20
circle_symbols_capacity = 9999
unique_yinyang_wuxing_symbols = 10  # 2 (yin/yang) * 5 (wuxing)

# "Book of Earth" - Square plate properties
square_side_cm = 10
square_symbols_capacity = 360
unique_bagua_symbols = 8

# Encoding target
num_chinese_chars_to_encode = 1000

# --- Step 2: Calculate characters per plate ---

# Symbols needed per character for circular plates (10^x >= 1000)
symbols_per_char_circle = math.ceil(math.log(num_chinese_chars_to_encode, unique_yinyang_wuxing_symbols))
# Characters per circle plate
chars_per_circle = math.floor(circle_symbols_capacity / symbols_per_char_circle)

# Symbols needed per character for squared plates (8^y >= 1000)
symbols_per_char_square = math.ceil(math.log(num_chinese_chars_to_encode, unique_bagua_symbols))
# Characters per square plate
chars_per_square = math.floor(square_symbols_capacity / symbols_per_char_square)


# --- Step 3: Solve the packing problem ---

# A circular plate requires a square of material with side equal to the diameter
circle_block_side_cm = circle_radius_cm * 2

# We discretize the material into a grid based on the smallest unit (10cm)
grid_unit_cm = square_side_cm
grid_w = material_width_cm // grid_unit_cm
grid_h = material_height_cm // grid_unit_cm

# Dimensions of the required material blocks in grid units
circle_block_grid_side = circle_block_side_cm // grid_unit_cm
square_block_grid_side = square_side_cm // grid_unit_cm

# To maximize characters, we prioritize the plate type with higher character density.
# Density (circle) = 3333 chars / (4*4 cells) = 208.3 chars/cell
# Density (square) = 90 chars / (1*1 cell) = 90 chars/cell
# Since circle density is higher, we maximize the number of circles first.

# Calculate the maximum number of circle blocks (4x4) that fit into the grid (14x11)
num_circles_along_w = grid_w // circle_block_grid_side
num_circles_along_h = grid_h // circle_block_grid_side
N_circles = num_circles_along_w * num_circles_along_h

# Calculate the number of square plates from the remaining area
total_cells = grid_w * grid_h
cells_used_by_circles = N_circles * (circle_block_grid_side ** 2)
remaining_cells = total_cells - cells_used_by_circles
M_squares = remaining_cells

# --- Step 4: Calculate total characters and format output ---

# K is the maximal number of characters
K_max_chars = N_circles * chars_per_circle + M_squares * chars_per_square

# Print the final answer and the equation used to derive it
print(f"To maximize the number of characters, Laozi should produce:")
print(f"- {N_circles} circular plates for the 'Book of Heaven'")
print(f"- {M_squares} squared plates for the 'Book of Earth'\n")

print("The calculation for the maximum number of characters (K) is:")
print(f"K = (Number of Circles * Chars/Circle) + (Number of Squares * Chars/Square)")
print(f"K = ({N_circles} * {chars_per_circle}) + ({M_squares} * {chars_per_square})")
print(f"K = {N_circles * chars_per_circle} + {M_squares * chars_per_square}")
print(f"K = {K_max_chars}\n")

print("The final answer in the N:M:K format is:")
print(f"{N_circles}:{M_squares}:{K_max_chars}")

<<<6:58:25218>>>