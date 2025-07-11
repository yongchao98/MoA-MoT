import math

# --- 1. Define Constants ---
# Tube dimensions
outer_side = 20  # cm
thickness = 4  # cm
tube_length = 100  # cm

# Ball dimensions
ball_radius = 2  # cm
ball_diameter = ball_radius * 2  # cm

# Values
value_whole_ball = 3
value_welded_ball = 2

# --- 2. Calculate Balls per Cross-Sectional Layer ---
# The cross-section is a frame. We can fit 4cm-diameter balls in it.
# Along the top and bottom outer edges (20cm):
num_balls_long_edge = math.floor(outer_side / ball_diameter)

# The side edges have a length of the inner dimension
inner_side = outer_side - 2 * thickness
# Along the two side edges (12cm):
num_balls_short_edge = math.floor(inner_side / ball_diameter)

# Total balls that can be placed in one layer of the cross-section
balls_per_layer = 2 * num_balls_long_edge + 2 * num_balls_short_edge

# --- 3. Evaluate Strategy 1: Making Whole Balls ---
# A whole ball requires a slice of material equal to its diameter
layers_for_whole_balls = math.floor(tube_length / ball_diameter)
total_whole_balls = balls_per_layer * layers_for_whole_balls
value_from_whole_balls = total_whole_balls * value_whole_ball

# --- 4. Evaluate Strategy 2: Making Half-Balls ---
# A half-ball requires a slice of material equal to its radius
layers_for_half_balls = math.floor(tube_length / ball_radius)
total_half_balls = balls_per_layer * layers_for_half_balls
# Two half-balls make one welded ball
total_welded_balls = math.floor(total_half_balls / 2)
value_from_welded_balls = total_welded_balls * value_welded_ball

# --- 5. Compare and Print the Best Strategy ---
print("To find the highest value, we analyze two strategies:")
print("\n--- Strategy 1: Making only whole balls ---")
print(f"The number of balls that fit in the cross-section is: 2 * floor({outer_side}/{ball_diameter}) + 2 * floor({inner_side}/{ball_diameter}) = {balls_per_layer}")
print(f"Number of {ball_diameter}cm layers in a {tube_length}cm tube: floor({tube_length} / {ball_diameter}) = {layers_for_whole_balls}")
print(f"Total whole balls = {balls_per_layer} balls/layer * {layers_for_whole_balls} layers = {total_whole_balls} balls")
print(f"Total value = {total_whole_balls} balls * {value_whole_ball} value/ball = {value_from_whole_balls}")

print("\n--- Strategy 2: Making only half-balls (welded) ---")
print(f"Number of {ball_radius}cm layers in a {tube_length}cm tube: floor({tube_length} / {ball_radius}) = {layers_for_half_balls}")
print(f"Total half-balls = {balls_per_layer} half-balls/layer * {layers_for_half_balls} layers = {total_half_balls} half-balls")
print(f"Total welded balls = floor({total_half_balls} / 2) = {total_welded_balls} balls")
print(f"Total value = {total_welded_balls} balls * {value_welded_ball} value/ball = {value_from_welded_balls}")

print("\n--- Conclusion ---")
if value_from_whole_balls > value_from_welded_balls:
    print(f"The highest possible value is {value_from_whole_balls} by making whole balls.")
    highest_value = value_from_whole_balls
else:
    print(f"The highest possible value is {value_from_welded_balls} by making welded balls.")
    highest_value = value_from_welded_balls
