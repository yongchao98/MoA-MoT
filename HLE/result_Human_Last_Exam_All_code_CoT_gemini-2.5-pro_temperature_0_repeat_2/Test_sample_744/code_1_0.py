import math

# Define the properties of the hollow tube and the balls
outer_side_cm = 20
length_m = 1
thickness_cm = 4
ball_radius_cm = 2

# Define the values
value_whole_ball = 3
value_welded_ball = 2

# --- Step 1: Determine the key dimensions ---
# Convert length to cm
length_cm = length_m * 100
# Calculate the ball diameter
ball_diameter_cm = ball_radius_cm * 2
# Calculate the inner side of the hollow tube
inner_side_cm = outer_side_cm - 2 * thickness_cm

print("Step 1: Calculating the number of balls that can be made from the material.")
print(f"The tube's material forms a frame with a cross-section of {outer_side_cm}x{outer_side_cm} cm and a hollow center of {inner_side_cm}x{inner_side_cm} cm.")
print(f"Each ball has a diameter of {ball_diameter_cm} cm.")

# --- Step 2: Calculate the number of balls ---
# Because the ball diameter (4cm) matches the tube's thickness (4cm), we can calculate the number of balls by dividing the volume.
# We can imagine the cross-section as a grid of 4x4 cm squares.
# Number of 4x4 squares in the outer 20x20 area
outer_grid_side = outer_side_cm // ball_diameter_cm
# Number of 4x4 squares in the inner 12x12 hollow area
inner_grid_side = inner_side_cm // ball_diameter_cm

# Number of 4x4 cm "slots" in the material's cross-section
slots_in_cross_section = outer_grid_side**2 - inner_grid_side**2
print(f"The number of balls that can fit in a single layer of the cross-section is ({outer_grid_side}^2 - {inner_grid_side}^2) = {slots_in_cross_section}.")

# Number of layers of balls along the tube's length
num_layers = length_cm // ball_diameter_cm
print(f"The number of layers along the {length_cm} cm length is {length_cm} / {ball_diameter_cm} = {num_layers}.")

# Total number of whole balls
total_whole_balls = slots_in_cross_section * num_layers
print(f"Total number of whole balls that can be manufactured = {slots_in_cross_section} * {num_layers} = {total_whole_balls}.")

# --- Step 3: Calculate the maximum value ---
# Since the value of a whole ball (3) is greater than a welded ball (2),
# and the material can be perfectly used to make whole balls, the optimal strategy is to only make whole balls.
max_value = total_whole_balls * value_whole_ball

print("\nStep 2: Calculating the highest possible value.")
print("To get the highest value, we should produce as many whole balls as possible, as they have a higher value (3) than welded balls (2).")
print("The final equation is the total number of whole balls multiplied by the value of a whole ball.")
print(f"{total_whole_balls} * {value_whole_ball} = {max_value}")
