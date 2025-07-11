import math

# Step 1: Define the dimensions
outer_side = 20  # cm
thickness = 4    # cm
length = 100     # cm (1m)
ball_radius = 2  # cm
ball_diameter = ball_radius * 2

# Step 2: Calculate how many balls can be cut from one cross-sectional layer
# The cross-section is a frame. We analyze its constituent rectangular parts.
# Top and bottom horizontal parts of the frame:
balls_in_horizontal_part = outer_side // ball_diameter

# Left and right vertical parts of the frame:
# The height of the vertical part available is the outer side minus the thickness occupied by the top and bottom balls.
vertical_space = outer_side - (2 * ball_diameter)
balls_in_vertical_part = vertical_space // ball_diameter

# Total balls per layer is the sum from all four sides of the frame.
balls_per_layer = (2 * balls_in_horizontal_part) + (2 * balls_in_vertical_part)

# Step 3: Calculate the number of layers along the tube's length
num_layers = length // ball_diameter

# Step 4: Calculate the total number of whole balls that can be manufactured
total_whole_balls = balls_per_layer * num_layers

# Step 5: Calculate the maximum value
# Value of a directly manufactured whole ball is 3.
# Value of a welded ball (from two half-balls) is 2.
# To maximize value, we should make as many high-value whole balls as possible.
# We assume the most straightforward interpretation: we cut the maximum number of whole balls,
# and the leftover material/shavings are discarded.
value_per_whole_ball = 3
max_value = total_whole_balls * value_per_whole_ball

# --- Output the results step-by-step ---
print(f"The tube has a wall thickness of {thickness} cm.")
print(f"Each ball has a diameter of {ball_diameter} cm.")
print("This means a ball can fit perfectly within the wall of the tube.")
print("\nCalculating balls per cross-section layer:")
print(f"  - Balls in top edge: floor({outer_side} / {ball_diameter}) = {balls_in_horizontal_part}")
print(f"  - Balls in bottom edge: floor({outer_side} / {ball_diameter}) = {balls_in_horizontal_part}")
print(f"  - Balls in each side edge: floor(({outer_side} - 2*{ball_diameter}) / {ball_diameter}) = {balls_in_vertical_part}")
print(f"Total balls per layer = {balls_in_horizontal_part} + {balls_in_horizontal_part} + {balls_in_vertical_part} + {balls_in_vertical_part} = {balls_per_layer}")

print(f"\nCalculating number of layers along the {length} cm length:")
print(f"Number of layers = floor({length} / {ball_diameter}) = {num_layers}")

print("\nCalculating total possible whole balls:")
print(f"Total whole balls = {balls_per_layer} (balls per layer) * {num_layers} (layers) = {total_whole_balls}")

print("\nCalculating the highest possible value:")
print(f"We prioritize the higher value item (whole balls at value 3).")
print(f"Highest Value = {total_whole_balls} balls * {value_per_whole_ball} value/ball")
print(f"Final Equation: {balls_per_layer} * {num_layers} * {value_per_whole_ball} = {max_value}")

<<<C>>>