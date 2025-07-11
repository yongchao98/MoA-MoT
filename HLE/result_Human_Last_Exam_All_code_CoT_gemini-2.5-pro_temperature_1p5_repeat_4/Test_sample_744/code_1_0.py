import math

# Step 1: Define the dimensions and values from the problem
# Tube dimensions
outer_side_cm = 20
length_cm = 100  # 1m = 100cm
thickness_cm = 4

# Ball dimensions
ball_radius_cm = 2

# Values
value_whole_ball = 3
value_welded_ball = 2

# Step 2: Calculate derived dimensions
inner_side_cm = outer_side_cm - 2 * thickness_cm
ball_diameter_cm = 2 * ball_radius_cm

print(f"The tube has a {outer_side_cm}x{outer_side_cm}cm outer cross-section and a {inner_side_cm}x{inner_side_cm}cm inner cross-section.")
print(f"The ball has a diameter of {ball_diameter_cm}cm.")
print("Since the wall thickness (4cm) equals the ball diameter (4cm), we can perfectly pack whole balls.")
print("-" * 30)

# Step 3: Calculate the number of balls that can be packed
# Calculate the cross-sectional area of the material
material_cross_section_area = outer_side_cm**2 - inner_side_cm**2

# Calculate the cross-sectional area of a ball's packing cube (4x4)
ball_packing_area = ball_diameter_cm**2

# Calculate how many rows of balls fit in the cross-section
num_rows = material_cross_section_area / ball_packing_area

# Calculate how many balls fit along the length of the tube
balls_per_row = math.floor(length_cm / ball_diameter_cm)

# Calculate the total number of whole balls
total_whole_balls = num_rows * balls_per_row

print(f"Number of rows that fit in the cross-section: {int(num_rows)}")
print(f"Number of balls that fit along the length: {balls_per_row}")
print(f"Total number of whole balls that can be made: {int(total_whole_balls)}")
print("-" * 30)


# Step 4: Calculate the maximum value
# Since whole balls are more valuable, we calculate the total value using only whole balls.
max_value = total_whole_balls * value_whole_ball

print("To maximize value, we should only make whole balls (value 3) instead of welded balls (value 2).")
print("\nFinal Calculation:")
print(f"{int(total_whole_balls)} balls * {value_whole_ball} value_per_ball = {int(max_value)}")
