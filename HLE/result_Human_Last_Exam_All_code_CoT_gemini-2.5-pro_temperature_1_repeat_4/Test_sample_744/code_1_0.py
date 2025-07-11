import math

# Step 1: Define the parameters of the hollow tube and the balls.
outer_side_cm = 20
thickness_cm = 4
length_cm = 100
ball_radius_cm = 2

value_whole_ball = 3
value_welded_ball = 2 # This value is for a ball made from 2 half-balls.

# Step 2: Calculate the cross-sectional area of the material.
inner_side_cm = outer_side_cm - 2 * thickness_cm
cross_section_area_material_cm2 = outer_side_cm**2 - inner_side_cm**2

# Step 3: Determine the optimal product and calculate how many can be made.
# As determined in the plan, making whole balls is the most valuable strategy.
ball_diameter_cm = ball_radius_cm * 2

# To cut a ball, we need a rod-like piece with a cross-section of ball_diameter x ball_diameter.
rod_cross_section_area_cm2 = ball_diameter_cm**2

# Calculate how many parallel rods can be cut from the tube's cross-section.
num_rods = cross_section_area_material_cm2 // rod_cross_section_area_cm2

# Calculate how many whole balls can be cut from the length of each rod.
num_whole_balls_per_rod = length_cm // ball_diameter_cm

# Calculate the total number of whole balls.
total_whole_balls = num_rods * num_whole_balls_per_rod

# Step 4: Calculate the final maximum value.
max_value = total_whole_balls * value_whole_ball

# Print the breakdown of the calculation as requested.
print("--- Calculation Breakdown ---")
print(f"1. Number of parallel 4x4cm rods that can be made from the tube: {int(num_rods)}")
print(f"2. Number of whole balls (4cm diameter) per 100cm rod: {int(num_whole_balls_per_rod)}")
print(f"3. Total number of whole balls produced: {int(num_rods)} rods * {int(num_whole_balls_per_rod)} balls/rod = {int(total_whole_balls)} balls")
print(f"4. Value of one directly manufactured whole ball: {value_whole_ball}")

print("\n--- Final Value Calculation ---")
print("Maximum Value = (Total Whole Balls) * (Value per Whole Ball)")
print("The final equation is:")
# The final part of the code prints each number in the final equation as requested.
print(f"{int(total_whole_balls)} * {value_whole_ball} = {int(max_value)}")
