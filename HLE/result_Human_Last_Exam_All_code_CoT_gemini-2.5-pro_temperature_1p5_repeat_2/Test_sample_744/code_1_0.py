import math

# Step 1: Define the dimensions
outer_dim = 20  # cm
length = 100  # cm (1m)
thickness = 4  # cm

ball_radius = 2  # cm
ball_diameter = ball_radius * 2

# Value of products
value_whole_ball = 3
value_welded_ball = 2

print("--- Problem Setup ---")
print(f"Hollow tube outer dimension: {outer_dim}x{outer_dim} cm")
print(f"Hollow tube length: {length} cm")
print(f"Hollow tube wall thickness: {thickness} cm")
print(f"Ball radius: {ball_radius} cm, Diameter: {ball_diameter} cm")
print("-" * 25)

# Step 2: Decompose the hollow tube into four rectangular prisms
# This decomposition avoids double-counting the corner material.
# We have two "top/bottom" prisms and two "side" prisms.

# Top and Bottom prisms
prism1_width = outer_dim
prism1_height = thickness
prism1_length = length

# Side prisms (material left after accounting for top and bottom)
inner_dim = outer_dim - 2 * thickness
prism2_width = inner_dim
prism2_height = thickness
prism2_length = length

print("--- Material Analysis ---")
print("The hollow tube's material can be seen as 4 rectangular prisms:")
print(f"2 prisms of size: {prism1_width} x {prism1_height} x {prism1_length} cm")
print(f"2 prisms of size: {prism2_width} x {prism2_height} x {prism2_length} cm")
print("-" * 25)


# Step 3: Calculate how many whole balls can be cut from each prism
# We are cutting balls with a 4cm diameter, which fit into 4x4x4 cm cubes.

# For the first pair of prisms (20x4x100)
balls_from_prism1 = (prism1_width // ball_diameter) * \
                    (prism1_height // ball_diameter) * \
                    (prism1_length // ball_diameter)
total_balls_from_prism1_pair = 2 * balls_from_prism1

# For the second pair of prisms (12x4x100)
balls_from_prism2 = (prism2_width // ball_diameter) * \
                    (prism2_height // ball_diameter) * \
                    (prism2_length // ball_diameter)
total_balls_from_prism2_pair = 2 * balls_from_prism2

total_whole_balls = total_balls_from_prism1_pair + total_balls_from_prism2_pair

print("--- Ball Manufacturing Calculation ---")
print(f"Balls from each of the two larger prisms: ({prism1_width}/{ball_diameter}) * ({prism1_height}/{ball_diameter}) * ({prism1_length}/{ball_diameter}) = {balls_from_prism1}")
print(f"Balls from each of the two smaller prisms: ({prism2_width}/{ball_diameter}) * ({prism2_height}/{ball_diameter}) * ({prism2_length}/{ball_diameter}) = {balls_from_prism2}")
print(f"Total whole balls: 2 * {balls_from_prism1} + 2 * {balls_from_prism2} = {total_whole_balls}")
print("-" * 25)

# Step 4: Analyze if any half-balls can be made
# Half-balls can be made from leftover material where a whole ball doesn't fit.
# A half-ball needs a piece of material at least 4x4x2 cm.
# We check for leftover material by using the modulo operator.
rem_prism1_width = prism1_width % ball_diameter
rem_prism1_height = prism1_height % ball_diameter
rem_prism1_length = prism1_length % ball_diameter

rem_prism2_width = prism2_width % ball_diameter
rem_prism2_height = prism2_height % ball_diameter
rem_prism2_length = prism2_length % ball_diameter

# Since all dimensions (20, 12, 4, 100) are perfectly divisible by the ball diameter (4),
# there is no leftover material to make half-balls.
welded_balls = 0

print("--- Half-Ball Analysis ---")
if rem_prism1_width == 0 and rem_prism1_height == 0 and rem_prism1_length == 0 and rem_prism2_width == 0 and rem_prism2_height == 0 and rem_prism2_length == 0:
    print("All material dimensions are perfectly divisible by the ball diameter.")
    print("There is no leftover material to create half-balls.")
    print(f"Number of welded balls: {welded_balls}")
print("-" * 25)


# Step 5: Calculate the maximum total value
total_value = (total_whole_balls * value_whole_ball) + (welded_balls * value_welded_ball)

print("--- Final Value Calculation ---")
print("The highest value is achieved by making the maximum number of whole balls, as they have a higher value.")
print(f"Final Equation: ({total_whole_balls} * {value_whole_ball}) + ({welded_balls} * {value_welded_ball}) = {total_value}")
print(f"The highest value I can make is: {total_value}")
