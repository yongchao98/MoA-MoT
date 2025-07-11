import math

# Define the properties of the hollow tube and the balls
outer_side_cm = 20
length_cm = 100  # 1m = 100cm
thickness_cm = 4
ball_radius_cm = 2

# Define the values for each type of ball
value_whole_ball = 3
value_welded_ball = 2

# Step 1: Calculate the total volume of the material
inner_side_cm = outer_side_cm - 2 * thickness_cm
outer_volume = outer_side_cm * outer_side_cm * length_cm
inner_volume = inner_side_cm * inner_side_cm * length_cm
material_volume = outer_volume - inner_volume

# Step 2: Determine the volume of material required to cut one ball
ball_diameter_cm = ball_radius_cm * 2
# Cutting a sphere requires a cube with side length equal to the sphere's diameter
volume_per_ball = ball_diameter_cm ** 3

# Step 3: Calculate the maximum number of balls that can be cut from the material
# As all material dimensions are multiples of the ball's diameter, we can use volume division
total_balls = material_volume / volume_per_ball

# Step 4: Choose the ball type with the maximum value
# Since value_whole_ball (3) > value_welded_ball (2), we always make whole balls
optimal_value_per_ball = value_whole_ball

# Step 5: Calculate the total maximum value
total_max_value = total_balls * optimal_value_per_ball

# Print the final equation and the result
print("To find the highest value, we first calculate the number of balls that can be made.")
print(f"The material volume is {int(material_volume)} cm³.")
print(f"Each ball requires a {int(ball_diameter_cm)}x{int(ball_diameter_cm)}x{int(ball_diameter_cm)} cube of material, which is {int(volume_per_ball)} cm³.")
print(f"Total number of balls = {int(material_volume)} / {int(volume_per_ball)} = {int(total_balls)}.")
print("To maximize value, we produce only whole balls, as they have a higher value (3) than welded balls (2).")
print("\nFinal Equation:")
print(f"{int(total_balls)} balls * {int(optimal_value_per_ball)} value/ball = {int(total_max_value)}")