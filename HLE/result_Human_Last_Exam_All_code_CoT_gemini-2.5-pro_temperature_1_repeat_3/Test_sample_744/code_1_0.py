import math

# Step 1: Define the dimensions and values from the problem statement.
outer_side_cm = 20
length_cm = 100 # 1m = 100cm
thickness_cm = 4
ball_radius_cm = 2
value_whole_ball = 3
value_welded_ball = 2

# Step 2: Calculate the volume of the material.
# The inner side is the outer side minus the thickness on both sides.
inner_side_cm = outer_side_cm - 2 * thickness_cm
# The volume is the cross-sectional area multiplied by the length.
cross_sectional_area = outer_side_cm**2 - inner_side_cm**2
material_volume_cm3 = cross_sectional_area * length_cm

print("Step 1: Calculating the volume of the raw material.")
print(f"The volume of the material is ({outer_side_cm}*{outer_side_cm} - {inner_side_cm}*{inner_side_cm}) * {length_cm} = {material_volume_cm3} cm^3.")
print("-" * 30)

# Step 3: Determine the number of whole balls that can be manufactured.
# The ball's diameter must fit within the material.
ball_diameter_cm = ball_radius_cm * 2
# Since the material thickness and ball diameter are both 4cm, we can cut whole balls.
# We can find the number of balls by seeing how many 4x4x4 "packing cubes" fit in the material volume.
packing_cube_volume = ball_diameter_cm ** 3
number_of_whole_balls = material_volume_cm3 / packing_cube_volume

print("Step 2: Determining the maximum number of balls.")
print(f"Each ball has a diameter of {ball_diameter_cm}cm.")
print("The number of balls is the total material volume divided by the volume of a 4x4x4 cube from which each ball is cut.")
print(f"Number of whole balls = {material_volume_cm3} / {packing_cube_volume} = {int(number_of_whole_balls)}")
print("-" * 30)

# Step 4: Calculate the maximum value.
# The optimal strategy is to make whole balls (value 3) instead of welded balls (value 2).
# Since the material can be entirely used for whole balls, there's no need to make any welded ones.
max_value = number_of_whole_balls * value_whole_ball

print("Step 3: Calculating the highest possible value.")
print("The highest value is achieved by making only whole balls, as they have a higher value (3) than welded balls (2).")
print("The final calculation is:")
print(f"{int(number_of_whole_balls)} * {value_whole_ball} = {int(max_value)}")

<<<C>>>