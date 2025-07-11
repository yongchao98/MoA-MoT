import math

# Step 1: Define the dimensions of the square hollow tube
outer_side_cm = 20
thickness_cm = 4
length_cm = 100 # 1m = 100cm

# Step 2: Define the properties of the balls
ball_radius_cm = 2
value_whole_ball = 3
value_welded_ball = 2

# Step 3: Calculate derived dimensions for packing
# The diameter of the ball determines the size of the 'cube' of material it occupies
ball_diameter_cm = ball_radius_cm * 2

# The inner side of the hollow tube's cross-section
inner_side_cm = outer_side_cm - 2 * thickness_cm

print(f"The tube has an outer side of {outer_side_cm}cm, an inner side of {inner_side_cm}cm, and a wall thickness of {thickness_cm}cm.")
print(f"Each ball has a diameter of {ball_diameter_cm}cm.")
print("Since the ball diameter is equal to the wall thickness, we can perfectly cut rows of whole balls.")
print("-" * 20)

# Step 4: Calculate how many balls can fit in a single cross-sectional slice
# Area of the material in the cross-section
material_cross_section_area = outer_side_cm**2 - inner_side_cm**2

# Area of the square bounding box required for one ball's cross-section
ball_box_area = ball_diameter_cm**2

# Number of balls per slice
balls_per_layer = material_cross_section_area // ball_box_area
print(f"The material's cross-sectional area is {outer_side_cm}^2 - {inner_side_cm}^2 = {material_cross_section_area} cm^2.")
print(f"Each ball requires a {ball_diameter_cm}x{ball_diameter_cm} area, which is {ball_box_area} cm^2.")
print(f"Number of balls per cross-sectional layer: {material_cross_section_area} / {ball_box_area} = {balls_per_layer}")
print("-" * 20)


# Step 5: Calculate how many layers of balls fit along the tube's length
num_layers = length_cm // ball_diameter_cm
print(f"The tube is {length_cm}cm long.")
print(f"Number of layers along the length: {length_cm}cm / {ball_diameter_cm}cm = {num_layers}")
print("-" * 20)

# Step 6: Calculate the total number of balls
total_balls = balls_per_layer * num_layers
print(f"Total number of balls that can be manufactured: {balls_per_layer} balls/layer * {num_layers} layers = {total_balls} balls.")
print("-" * 20)

# Step 7: Calculate the maximum possible value
# To maximize value, we should produce whole balls (value 3) instead of welded balls (value 2).
max_value = total_balls * value_whole_ball

print("To get the highest value, we make all 400 balls as whole balls.")
print(f"The highest possible value is:")
print(f"{total_balls} balls * {value_whole_ball} (value per whole ball) = {max_value}")

<<<C>>>