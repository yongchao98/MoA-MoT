import math

# --- Given Parameters ---
# Tube dimensions
outer_dim_cm = 20
thickness_cm = 4
length_cm = 100  # 1m

# Ball dimensions
ball_radius_cm = 2

# Product values
value_whole_ball = 3
value_welded_ball = 2

# --- Step-by-Step Calculation ---

# Step 1: Determine the bounding box for cutting a ball
ball_diameter_cm = ball_radius_cm * 2

# Step 2: Calculate the dimensions of the inner hollow space
inner_dim_cm = outer_dim_cm - 2 * thickness_cm

# Step 3: Calculate how many balls can be cut from a single cross-sectional slice (4cm thick)
# We find how many 4x4 squares fit in the outer area and subtract those in the hollow inner area.
n_outer = outer_dim_cm // ball_diameter_cm
n_inner = inner_dim_cm // ball_diameter_cm

balls_in_outer_area = n_outer ** 2
balls_in_inner_area = n_inner ** 2
balls_per_slice = balls_in_outer_area - balls_in_inner_area

# Step 4: Calculate how many 4cm-thick slices we can get from the tube's length
num_slices = length_cm // ball_diameter_cm

# Step 5: Calculate the total number of balls that can be manufactured
# Since a whole ball (value 3) is more valuable than a welded ball (value 2) and
# there's no geometric reason to make half-balls, we will only make whole balls.
num_whole_balls = balls_per_slice * num_slices
num_welded_balls = 0

# Step 6: Calculate the total maximum value
total_value = (num_whole_balls * value_whole_ball) + (num_welded_balls * value_welded_ball)

# --- Output the result ---
print("This problem is about packing, not just dividing volume.")
print(f"A ball with radius {ball_radius_cm}cm has a diameter of {ball_diameter_cm}cm.")
print("We can determine how many balls can be cut based on these dimensions.\n")
print("1. Balls per cross-section slice:")
print(f"   - The 20x20cm outer section can fit {n_outer}x{n_outer} = {balls_in_outer_area} balls.")
print(f"   - The {inner_dim_cm}x{inner_dim_cm}cm inner hollow section removes {n_inner}x{n_inner} = {balls_in_inner_area} balls.")
print(f"   - Net balls per 4cm-thick slice = {balls_in_outer_area} - {balls_in_inner_area} = {balls_per_slice} balls.\n")
print("2. Number of slices along the 100cm length:")
print(f"   - Number of 4cm slices = {length_cm}cm / {ball_diameter_cm}cm = {num_slices} slices.\n")
print("3. Total balls possible:")
print(f"   - Total balls = {balls_per_slice} balls/slice * {num_slices} slices = {num_whole_balls} balls.\n")
print("4. Calculating maximum value:")
print(f"   - Since whole balls (value 3) are more valuable than welded balls (value 2), we make only whole balls.")
print("   - Final Equation:")
print(f"   - Total Value = {num_whole_balls} (whole balls) * {value_whole_ball} (value) + {num_welded_balls} (welded balls) * {value_welded_ball} (value) = {total_value}")

print(f"\nThe highest value I can make from my material is {total_value}.")