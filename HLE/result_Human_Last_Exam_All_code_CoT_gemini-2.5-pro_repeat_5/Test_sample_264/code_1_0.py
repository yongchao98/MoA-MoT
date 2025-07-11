import sys

# Step 1: Define the properties of the cubes and the shift.
# The side length of the original cubes.
cube_side_length = 1.0
# The shift along each axis.
shift_distance = 0.5

# The first cube's boundaries are [0, 1] on each axis.
cube1_min = 0.0
cube1_max = cube_side_length

# The second cube's boundaries are shifted by 0.5 on each axis.
cube2_min = shift_distance
cube2_max = cube_side_length + shift_distance

# Step 2: Find the boundaries of the intersection.
# The intersection is a smaller cube defined by the overlapping region on each axis.
# For each axis, the overlap is from max(cube1_min, cube2_min) to min(cube1_max, cube2_max).
intersection_min = max(cube1_min, cube2_min)
intersection_max = min(cube1_max, cube2_max)

# Step 3: Calculate the side length of the intersection cube.
intersection_side_length = intersection_max - intersection_min

# Step 4: Calculate the volume of the intersection.
volume = intersection_side_length ** 3

# Print the explanation and the final calculation.
print(f"The side length of the intersection cube is the length of the overlapping region on one axis.")
print(f"Side Length = {intersection_max} - {intersection_min} = {intersection_side_length}")
print(f"\nThe volume of the intersection is the cube of its side length.")
# The final equation with each number.
print(f"Volume = {intersection_side_length} * {intersection_side_length} * {intersection_side_length} = {volume}")

# The user is running the script, so we write the answer to the file 'answer.txt'
# in the format they requested.
try:
    with open('answer.txt', 'w') as f:
        f.write(f'<<<{volume}>>>')
except (IOError, OSError) as e:
    # Exit gracefully if we can't write the file.
    # We still printed the result to the console.
    sys.exit()