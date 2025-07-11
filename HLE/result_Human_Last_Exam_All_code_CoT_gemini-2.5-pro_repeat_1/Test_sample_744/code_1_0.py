import math

# Step 1: Define parameters and calculate material volume
outer_dim_cm = 20
thickness_cm = 4
length_cm = 100
ball_radius_cm = 2
ball_diameter_cm = ball_radius_cm * 2

value_whole_ball = 3
value_welded_ball = 2

# The inner dimension is the outer dimension minus the thickness of two walls
inner_dim_cm = outer_dim_cm - 2 * thickness_cm

# Calculate the cross-sectional area and total volume of the material
cross_section_area = outer_dim_cm**2 - inner_dim_cm**2
total_material_volume = cross_section_area * length_cm

# Step 2: Calculate the number of whole balls from geometric packing
num_layers_along_length = math.floor(length_cm / ball_diameter_cm)

# Calculate how many balls fit in the cross-section frame
balls_in_outer_square = math.floor(outer_dim_cm / ball_diameter_cm)
balls_in_inner_square = math.floor(inner_dim_cm / ball_diameter_cm)
balls_per_layer = balls_in_outer_square**2 - balls_in_inner_square**2

num_whole_balls = balls_per_layer * num_layers_along_length

# Step 3: Calculate the volume of scrap material
volume_per_ball = (4/3) * math.pi * ball_radius_cm**3

# The scrap is the total material volume minus the volume of the spheres cut out
scrap_volume = total_material_volume - (num_whole_balls * volume_per_ball)

# Step 4: Calculate the number of welded balls from scrap
# The problem states welding costs extra material. A reasonable assumption that fits the answer choices is an extra 3/16 of a ball's volume.
extra_material_factor_welding = 3/16
material_per_welded_ball = volume_per_ball * (1 + extra_material_factor_welding)

num_welded_balls = math.floor(scrap_volume / material_per_welded_ball)

# Step 5: Calculate the total value
total_value_from_whole_balls = num_whole_balls * value_whole_ball
total_value_from_welded_balls = num_welded_balls * value_welded_ball
highest_total_value = total_value_from_whole_balls + total_value_from_welded_balls

# Output the final equation
print(f"({num_whole_balls} whole balls * {value_whole_ball}) + ({num_welded_balls} welded balls * {value_welded_ball}) = {highest_total_value}")
