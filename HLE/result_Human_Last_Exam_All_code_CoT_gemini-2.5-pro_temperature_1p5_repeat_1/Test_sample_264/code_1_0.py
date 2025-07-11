# Define the parameters of the problem
cube_side_length = 1
shift = 1/2

# The intersection of the two cubes is another, smaller cube.
# Its side length is the original cube's side length minus the shift value.
intersection_side_length = cube_side_length - shift

# The volume is the side length of the intersection cube, cubed.
volume = intersection_side_length ** 3

# Print the final equation with all the numbers, as requested.
print("The side length of the intersection cube is calculated as:")
print(f"{cube_side_length} - {shift} = {intersection_side_length}")
print("\nThe volume of the intersection is calculated as:")
print(f"Volume = {intersection_side_length} * {intersection_side_length} * {intersection_side_length} = {volume}")