import math

# Step 1: Define the parameters
outer_side = 20  # cm
thickness = 4  # cm
length = 100  # cm (1m)
ball_radius = 2  # cm
value_whole_ball = 3
value_welded_ball = 2

# Step 2: Calculate the cross-sectional area of the material
inner_side = outer_side - 2 * thickness
outer_area = outer_side ** 2
inner_area = inner_side ** 2
material_cross_section_area = outer_area - inner_area

print(f"The cross-sectional area of the material is {material_cross_section_area} cm^2.")

# Step 3: Determine ball packing
ball_diameter = 2 * ball_radius
print(f"The ball diameter is {ball_diameter} cm.")

# The wall thickness is 4 cm, same as the ball diameter.
# We can find how many ball 'slots' (4x4 cm) fit in the cross-section.
num_stacks_in_cross_section = material_cross_section_area / (ball_diameter ** 2)
print(f"Number of stacks that can fit in the cross-section: {int(num_stacks_in_cross_section)}")

# Step 4: Calculate balls per stack
num_balls_per_stack = length / ball_diameter
print(f"Number of balls per stack along the length: {int(num_balls_per_stack)}")

# Step 5: Calculate total number of whole balls
total_whole_balls = num_stacks_in_cross_section * num_balls_per_stack
print(f"Total number of whole balls that can be made: {int(total_whole_balls)}")

# Step 6: Determine the optimal strategy and calculate the maximum value
# Since a whole ball (value 3) is more valuable than a welded ball (value 2),
# and the dimensions allow for manufacturing with no leftover material,
# the best strategy is to make only whole balls.
max_value = total_whole_balls * value_whole_ball
print("\nFinal Calculation:")
print(f"The highest value is obtained by making only whole balls, as there is no leftover material.")
print(f"Total Value = {int(total_whole_balls)} (balls) * {value_whole_ball} (value per ball) = {int(max_value)}")

<<<C>>>