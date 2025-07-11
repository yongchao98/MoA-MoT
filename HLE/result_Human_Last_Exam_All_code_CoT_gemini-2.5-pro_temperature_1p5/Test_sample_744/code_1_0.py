import math

# Step 1: Define the dimensions
tube_outer_side = 20  # cm
tube_length = 100  # cm, since 1m = 100cm
tube_thickness = 4  # cm
ball_radius = 2  # cm
ball_diameter = ball_radius * 2

# Step 2: Define values
value_whole_ball = 3
value_welded_ball = 2

# Step 3: Calculate the number of balls from the corner pieces
# The cross-section of each corner piece is thickness x thickness (4x4 cm)
# The length is the tube length (100 cm)
# We find how many balls (diameter 4cm) fit along each dimension
num_balls_per_corner_length = math.floor(tube_length / ball_diameter)
num_balls_per_corner_width = math.floor(tube_thickness / ball_diameter)
num_balls_per_corner_height = math.floor(tube_thickness / ball_diameter)

# Total balls in one corner piece
balls_in_one_corner = num_balls_per_corner_length * num_balls_per_corner_width * num_balls_per_corner_height
# Total balls from all 4 corners
total_balls_from_corners = 4 * balls_in_one_corner

print(f"Calculation for corner pieces:")
print(f"Each of the 4 corner pieces is {tube_thickness}x{tube_thickness}x{tube_length} cm.")
print(f"Number of balls that fit in one corner piece = floor({tube_length}/{ball_diameter}) * floor({tube_thickness}/{ball_diameter}) * floor({tube_thickness}/{ball_diameter}) = {num_balls_per_corner_length} * {num_balls_per_corner_width} * {num_balls_per_corner_height} = {balls_in_one_corner}")
print(f"Total balls from 4 corner pieces = 4 * {balls_in_one_corner} = {total_balls_from_corners}\n")


# Step 4: Calculate the number of balls from the side pieces
# The side plank width is the tube's inner side length
inner_side_length = tube_outer_side - 2 * tube_thickness
# The side plank cross-section is inner_side_length x thickness (12x4 cm)
num_balls_per_side_length = math.floor(tube_length / ball_diameter)
num_balls_per_side_width = math.floor(inner_side_length / ball_diameter)
num_balls_per_side_height = math.floor(tube_thickness / ball_diameter)

# Total balls in one side piece
balls_in_one_side = num_balls_per_side_length * num_balls_per_side_width * num_balls_per_side_height
# Total balls from all 4 sides
total_balls_from_sides = 4 * balls_in_one_side

print(f"Calculation for side pieces:")
print(f"The inner side length is {tube_outer_side} - 2*{tube_thickness} = {inner_side_length} cm.")
print(f"Each of the 4 side pieces is {inner_side_length}x{tube_thickness}x{tube_length} cm.")
print(f"Number of balls that fit in one side piece = floor({tube_length}/{ball_diameter}) * floor({inner_side_length}/{ball_diameter}) * floor({tube_thickness}/{ball_diameter}) = {num_balls_per_side_length} * {num_balls_per_side_width} * {num_balls_per_side_height} = {balls_in_one_side}")
print(f"Total balls from 4 side pieces = 4 * {balls_in_one_side} = {total_balls_from_sides}\n")

# Step 5: Calculate total balls and maximum value
total_whole_balls = total_balls_from_corners + total_balls_from_sides
max_value = total_whole_balls * value_whole_ball

print("Final Calculation:")
print("The material dimensions are perfect multiples of the ball diameter.")
print("This means we can only produce whole balls, as there are no leftover sections for half-balls.")
print(f"Total whole balls = Balls from corners + Balls from sides = {total_balls_from_corners} + {total_balls_from_sides} = {total_whole_balls}")
print(f"The highest value is obtained by making only whole balls.")
print(f"Highest Value = Total Whole Balls * Value per Whole Ball")
print(f"Highest Value = {total_whole_balls} * {value_whole_ball} = {max_value}")
