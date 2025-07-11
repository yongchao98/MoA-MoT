# Define the dimensions of the material and the product
outer_side_cm = 20
length_cm = 100
thickness_cm = 4
ball_radius_cm = 2

# Calculate the diameter of the ball, which determines the size of the cube needed
ball_diameter_cm = ball_radius_cm * 2

# Calculate the dimensions of the inner hollow space
inner_side_cm = outer_side_cm - 2 * thickness_cm

# Calculate how many 4x4x4 cubes can be made from the material
# First, find the number of cubes in the outer solid prism
cubes_in_outer_prism = (outer_side_cm // ball_diameter_cm) * (outer_side_cm // ball_diameter_cm) * (length_cm // ball_diameter_cm)

# Second, find the number of cubes in the inner hollow prism
cubes_in_inner_prism = (inner_side_cm // ball_diameter_cm) * (inner_side_cm // ball_diameter_cm) * (length_cm // ball_diameter_cm)

# The total number of balls is the difference
total_whole_balls = cubes_in_outer_prism - cubes_in_inner_prism

# Define the values for each type of ball
value_whole_ball = 3
value_welded_ball = 2

# To maximize the value, we should make as many whole balls as possible,
# as they have a higher value (3) than welded balls (2) for the same amount of material.
# Since the material dimensions are perfectly divisible by the ball diameter,
# we can make 400 whole balls with no leftover material to be used for half-balls.
max_value = total_whole_balls * value_whole_ball

# Print the final equation
print(f"The highest value is obtained by making {total_whole_balls} whole balls.")
print(f"The calculation is:")
print(f"{total_whole_balls} * {value_whole_ball} = {max_value}")
