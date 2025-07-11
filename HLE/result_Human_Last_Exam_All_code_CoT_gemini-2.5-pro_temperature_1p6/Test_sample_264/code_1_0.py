# Side length of the original cubes
side_length = 1.0

# The shift value along each axis. "Shifted by 1/2 along the main diagonal"
# is interpreted as a shift vector of (1/2, 1/2, 1/2).
shift_per_axis = 0.5

# Calculate the side length of the intersection region.
# On any given axis, the first cube is at [0, 1] and the second is at [0.5, 1.5].
# Their intersection is the interval [0.5, 1], which has a length of 0.5.
intersection_side = side_length - shift_per_axis

# The volume is the side length of the intersection cube, cubed.
volume = intersection_side ** 3

# Print the final equation showing the multiplication of the side lengths.
print(f"The side length of the intersection cube is {intersection_side}.")
print(f"The volume calculation is:")
print(f"{intersection_side} * {intersection_side} * {intersection_side} = {volume}")