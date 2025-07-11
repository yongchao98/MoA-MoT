import math

# Step 1: Define the dimensions of the shapes based on the problem's constraints.
r_yellow = 1
# From the geometric analysis, the radius of a white circle is 4.
r_white = 4
d_white = 2 * r_white

# Step 2: Determine the overall dimensions of the target.
# Green shapes are squares with side length equal to the diameter of a white circle.
w_green = d_white
h_green = d_white

# The total width is determined by the rows, which contain white circles and green squares.
# E.g., Middle Row: 1 green square + 4 white circles = 8 + 4*8 = 40
num_white_circles_per_row = 4
num_green_squares_per_row = 1
total_width = num_white_circles_per_row * d_white + num_green_squares_per_row * w_green

# The total height is determined by the 3 rows of white circles.
num_rows = 3
total_height = num_rows * d_white

# Calculate the total area of the target.
total_area = total_width * total_height
print(f"The total width of the target is {total_width} cm.")
print(f"The total height of the target is {total_height} cm.")
print(f"The total area of the target is {total_width} * {total_height} = {total_area} cm^2.\n")


# Step 3: Find the total area of the yellow circles.
# Based on our analysis, there are 12 yellow circles on the target.
# 4 on the top edge, 4 on the bottom edge, 2 on the left edge, 2 on the right edge.
num_yellow_circles = 12
area_one_yellow_circle = math.pi * (r_yellow ** 2)
total_yellow_area = num_yellow_circles * area_one_yellow_circle

print(f"There are {num_yellow_circles} yellow circles.")
print(f"The area of one yellow circle is pi * {r_yellow}^2 = {area_one_yellow_circle:.4f} cm^2.")
print(f"The total area of all yellow circles is {num_yellow_circles} * pi = {total_yellow_area:.4f} cm^2.\n")


# Step 4: Calculate the expected number of hits in 10000 shots.
num_shots = 10000
probability_hit_yellow = total_yellow_area / total_area
expected_hits = num_shots * probability_hit_yellow

print("The probability of hitting a yellow circle is (Total Yellow Area) / (Total Target Area).")
print(f"P = ({num_yellow_circles} * pi) / {total_area} = pi / {int(total_area/num_yellow_circles)}")
print(f"Expected hits in {num_shots} shots = {num_shots} * (({num_yellow_circles} * pi) / {total_area})")
print(f"Expected hits = {expected_hits:.4f}")

# Final equation representation
print("\nFinal Equation:")
print(f"Expected Hits = {num_shots} * ({num_yellow_circles} * \u03C0 * {r_yellow}\u00b2) / ({total_width} * {total_height}) = {expected_hits:.4f}")

# The final answer in the required format
final_answer_value = 125 * math.pi
print(f"\nFinal Answer Value: {final_answer_value}")
<<<392.69908169872417>>>