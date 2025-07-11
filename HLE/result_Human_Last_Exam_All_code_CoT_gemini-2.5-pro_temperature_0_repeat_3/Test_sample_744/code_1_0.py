import math

# --- Problem Parameters ---

# Tube dimensions in cm
outer_side = 20
thickness = 4
length = 100  # 1m = 100cm

# Ball dimensions in cm
ball_radius = 2

# Values
value_whole_ball = 3
value_welded_ball = 2

# --- Calculations ---

# 1. Calculate the space required for one ball
# A ball with a 2cm radius has a 4cm diameter.
# We assume it's cut from a 4x4x4 cm cube of material.
ball_diameter = ball_radius * 2

# 2. Calculate how many balls can be packed into the material
# First, calculate the dimensions of the inner hollow space
inner_side = outer_side - 2 * thickness

# Calculate the cross-sectional area of the material
material_cross_section_area = outer_side**2 - inner_side**2

# Calculate the cross-sectional area of the bounding box for one ball
ball_bounding_box_area = ball_diameter**2

# The number of balls that can fit in a single 4cm thick slice of the tube
# This works because the tube's thickness (4cm) is the same as the ball's diameter (4cm)
num_balls_per_layer = material_cross_section_area / ball_bounding_box_area

# Calculate the number of layers of balls that fit along the tube's length
num_layers = length / ball_diameter

# Total number of locations where a ball can be manufactured
total_ball_locations = num_balls_per_layer * num_layers

# 3. Determine the optimal strategy and calculate the maximum value
# At each location, we can make a whole ball (value 3) or a welded ball (value 2).
# To maximize value, we should always choose the option with the higher value.
# Since 3 > 2, we will only make whole balls.
max_value = total_ball_locations * value_whole_ball

# --- Output the result ---
# The problem asks to output each number in the final equation.
print("The highest value can be calculated by multiplying the total number of whole balls by the value of a single whole ball.")
print("Final Equation:")
print(f"{int(total_ball_locations)} * {value_whole_ball} = {int(max_value)}")
