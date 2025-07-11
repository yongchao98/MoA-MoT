import math

# Step 1: Define the parameters from the problem description
tube_outer_width = 20  # cm
tube_length = 100  # cm, since 1m = 100cm
tube_thickness = 4  # cm
ball_radius = 2  # cm

# Define the value of the products
whole_ball_value = 3
welded_ball_value = 2

# Step 2: Calculate derived dimensions
ball_diameter = ball_radius * 2
# The inner width is the outer width minus the thickness of the two walls
tube_inner_width = tube_outer_width - 2 * tube_thickness

print("--- Calculating Dimensions ---")
print(f"The material is a square hollow tube with an outer width of {tube_outer_width} cm, a thickness of {tube_thickness} cm, and a length of {tube_length} cm.")
print(f"The balls to be made have a radius of {ball_radius} cm, which means a diameter of {ball_diameter} cm.")
print(f"The tube's inner width is {tube_outer_width} - 2 * {tube_thickness} = {tube_inner_width} cm.")
print("\n")

# Step 3: Calculate how many balls can be packed into the cross-section
# We can visualize the cross-section as being made of 4 side rectangles and 4 corner squares.
# Side rectangle dimensions are (tube_inner_width) x (tube_thickness) -> 12cm x 4cm
# Corner square dimensions are (tube_thickness) x (tube_thickness) -> 4cm x 4cm

# Number of balls in the cross-section of one side piece
balls_in_side_cross_section = (tube_inner_width // ball_diameter) * (tube_thickness // ball_diameter)

# Number of balls in the cross-section of one corner piece
balls_in_corner_cross_section = (tube_thickness // ball_diameter) * (tube_thickness // ball_diameter)

# Total balls per cross-sectional layer
total_balls_per_layer = 4 * balls_in_side_cross_section + 4 * balls_in_corner_cross_section

print("--- Calculating Packing ---")
print("We assume a simple cubic packing arrangement since the ball diameter (4cm) fits perfectly into the wall thickness (4cm).")
print(f"Number of balls that fit in one side-rectangle's cross-section (12x4 cm): ({tube_inner_width}//{ball_diameter}) * ({tube_thickness}//{ball_diameter}) = {balls_in_side_cross_section}")
print(f"Number of balls that fit in one corner-square's cross-section (4x4 cm): ({tube_thickness}//{ball_diameter}) * ({tube_thickness}//{ball_diameter}) = {balls_in_corner_cross_section}")
print(f"Total balls per layer = 4 * {balls_in_side_cross_section} (sides) + 4 * {balls_in_corner_cross_section} (corners) = {total_balls_per_layer} balls.")
print("\n")

# Step 4: Calculate how many layers of balls fit along the tube's length
num_layers = tube_length // ball_diameter
print(f"Number of layers along the tube's length ({tube_length} cm): {tube_length} // {ball_diameter} = {num_layers} layers.")
print("\n")

# Step 5: Calculate total number of manufacturable balls
# Since we prioritize making more valuable whole balls, we calculate the total number of whole balls.
# The dimensions allow for perfect packing, so there is no leftover material to make half-balls from.
total_whole_balls = total_balls_per_layer * num_layers

print("--- Calculating Total Balls ---")
print(f"Total number of whole balls = (Balls per Layer) * (Number of Layers)")
print(f"Total number of whole balls = {total_balls_per_layer} * {num_layers} = {total_whole_balls}")
print("\n")

# Step 6: Calculate the maximum value
# The highest value is achieved by making only whole balls.
max_value = total_whole_balls * whole_ball_value

print("--- Calculating Maximum Value ---")
print("To get the highest value, we make as many whole balls as possible.")
print("Each whole ball has a value of 3.")
print("The information about welded balls (value 2) is not needed as the tube's dimensions are perfect multiples of the ball's diameter, leaving no geometric scraps to be turned into half-balls.")
print("\nFinal equation:")
print(f"{total_whole_balls} * {whole_ball_value} = {max_value}")
