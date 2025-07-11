# Define the properties of the material and artifacts
sheet_width = 140
sheet_height = 110
circle_diameter = 40
circle_value = 9999
square_side = 10
square_value = 360

# --- Step 1: Maximize the number of high-value circular plates ---
# Calculate how many circles fit in a grid layout on the sheet
num_circles_along_width = sheet_width // circle_diameter
num_circles_along_height = sheet_height // circle_diameter
num_circles = num_circles_along_width * num_circles_along_height

print(f"To maximize characters, we first fit the maximum number of circular plates.")
print(f"Number of circular plates (N): {num_circles}")

# --- Step 2: Calculate space for squared plates from the remainder ---
# Calculate the dimensions of the rectangular block used by the circles
circles_used_width = num_circles_along_width * circle_diameter
circles_used_height = num_circles_along_height * circle_diameter

# The remaining space is an L-shape, which we split into two rectangles
# to calculate the number of squares that can fit.
# Rectangle 1:
rem_1_w = sheet_width
rem_1_h = sheet_height - circles_used_height
sq_in_rem_1 = (rem_1_w // square_side) * (rem_1_h // square_side)

# Rectangle 2:
rem_2_w = sheet_width - circles_used_width
rem_2_h = circles_used_height
sq_in_rem_2 = (rem_2_w // square_side) * (rem_2_h // square_side)

num_squares = sq_in_rem_1 + sq_in_rem_2

print(f"Number of squared plates (M) from remaining material: {num_squares}")

# --- Step 3: Calculate the maximum total number of characters ---
max_characters = (num_circles * circle_value) + (num_squares * square_value)

print("\nThe final calculation for the maximum number of characters is:")
print(f"{num_circles} * {circle_value} + {num_squares} * {square_value} = {max_characters}")

# --- Final Answer ---
# Present the final answer in the format N:M:K
final_answer_format = f"{num_circles}:{num_squares}:{max_characters}"
print(f"\nThe optimal configuration is N:M:K = {final_answer_format}")

# Use the required output format for the final answer
print(f"\n<<<{final_answer_format}>>>")