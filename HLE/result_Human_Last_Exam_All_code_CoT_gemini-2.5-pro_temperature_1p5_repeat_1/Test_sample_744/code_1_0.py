import math

# Define constants based on the problem statement
tube_outer_side = 20  # cm
tube_length = 100  # cm (1m)
tube_thickness = 4  # cm
ball_radius = 2  # cm

value_whole_ball = 3
value_welded_ball = 2

# Derived dimensions needed for calculation
ball_diameter = ball_radius * 2
ball_height = ball_diameter  # for a whole ball
half_ball_height = ball_radius

print("--- Step 1: Analyze Strategy of Making Whole Balls ---")

# Calculate how many balls fit in the cross-section
# A 4cm diameter ball has a 4x4cm footprint.
# We find how many 4x4cm squares fit in the tube's cross-section.
outer_grid_side = tube_outer_side / ball_diameter
inner_tube_side = tube_outer_side - 2 * tube_thickness
inner_grid_side = inner_tube_side / ball_diameter

balls_per_layer = int(outer_grid_side**2 - inner_grid_side**2)
print(f"Number of balls that can fit in the cross-section per layer: {balls_per_layer}")

# Calculate how many layers of whole balls fit along the length
num_layers_whole_balls = int(tube_length / ball_height)
print(f"Number of layers of whole balls along the 100cm length: {num_layers_whole_balls}")

# Calculate total whole balls and value
total_whole_balls = balls_per_layer * num_layers_whole_balls
value_from_whole_balls = total_whole_balls * value_whole_ball
print(f"Total whole balls: {balls_per_layer} * {num_layers_whole_balls} = {total_whole_balls}")
print(f"Equation for total value from whole balls: {total_whole_balls} * {value_whole_ball} = {value_from_whole_balls}")

print("\n--- Step 2: Analyze Strategy of Making Welded Balls from Half-Balls ---")

# The cross-section footprint of a half-ball is the same as a whole ball
half_balls_per_layer = balls_per_layer
print(f"Number of half-balls that can fit in the cross-section per layer: {half_balls_per_layer}")

# Calculate how many layers of half-balls fit along the length
num_layers_half_balls = int(tube_length / half_ball_height)
print(f"Number of layers of half-balls along the 100cm length: {num_layers_half_balls}")

# Calculate total half-balls and the number of resulting welded balls
total_half_balls = half_balls_per_layer * num_layers_half_balls
total_welded_balls = int(total_half_balls / 2)
print(f"Total half-balls: {half_balls_per_layer} * {num_layers_half_balls} = {total_half_balls}")
print(f"Total welded balls: {total_half_balls} / 2 = {total_welded_balls}")

# Calculate the total value from welded balls
value_from_welded_balls = total_welded_balls * value_welded_ball
print(f"Equation for total value from welded balls: {total_welded_balls} * {value_welded_ball} = {value_from_welded_balls}")


print("\n--- Step 3: Conclusion ---")
# Compare the values from the two strategies to find the maximum
if value_from_whole_balls > value_from_welded_balls:
    highest_value = value_from_whole_balls
    print(f"Comparing the two strategies ({value_from_whole_balls} vs {value_from_welded_balls}), making whole balls yields a higher value.")
else:
    highest_value = value_from_welded_balls
    print(f"Comparing the two strategies ({value_from_whole_balls} vs {value_from_welded_balls}), making welded balls yields a higher value.")

print(f"\nThe highest value I can make from my material is: {highest_value}")