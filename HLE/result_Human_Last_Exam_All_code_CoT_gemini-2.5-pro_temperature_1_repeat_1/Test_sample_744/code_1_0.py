# Define the properties of the material and the product
outer_side_cm = 20
thickness_cm = 4
length_cm = 100  # 1 meter

ball_radius_cm = 2

value_whole_ball = 3
value_welded_ball = 2 # This is lower, so we prioritize whole balls

# --- Step 1: Calculate the usable cross-sectional area of the tube ---

# The inner side is the outer side minus the thickness on both sides
inner_side_cm = outer_side_cm - 2 * thickness_cm

# The area is the outer square's area minus the inner hollow square's area
material_cross_section_area_cm2 = outer_side_cm**2 - inner_side_cm**2

# --- Step 2: Calculate how many balls can be cut ---

# The ball's diameter determines its cutting profile
ball_diameter_cm = ball_radius_cm * 2

# The cutting profile of a ball is a square with side equal to its diameter
ball_profile_area_cm2 = ball_diameter_cm**2

# The number of parallel columns of balls we can cut
# The material thickness (4cm) perfectly fits the ball diameter (4cm),
# so we can use the entire cross-sectional area.
num_columns = material_cross_section_area_cm2 // ball_profile_area_cm2

# The number of balls that fit along the 1m (100cm) length
num_balls_per_column = length_cm // ball_diameter_cm

# Total number of whole balls that can be manufactured
total_balls = num_columns * num_balls_per_column

# --- Step 3: Calculate the highest possible value ---

# The highest value is achieved by making whole balls, which are more valuable.
max_value = total_balls * value_whole_ball

# --- Step 4: Print the final result in an equation format ---

print(f"The number of parallel columns of balls that can be cut is {num_columns}.")
print(f"The number of balls per column is {num_balls_per_column}.")
print(f"Total number of whole balls that can be manufactured: {num_columns} * {num_balls_per_column} = {total_balls}")
print(f"The highest possible value is calculated as:")
print(f"{total_balls} (total balls) * {value_whole_ball} (value per whole ball) = {max_value}")
