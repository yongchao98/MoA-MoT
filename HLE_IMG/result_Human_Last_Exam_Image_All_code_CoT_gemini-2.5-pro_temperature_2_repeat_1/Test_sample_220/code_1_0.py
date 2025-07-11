import math

# A face-centered cubic (FCC) crystal structure, when viewed along the [110] direction,
# projects its atomic positions into a 2D pattern known as a centered rectangle.
# This script calculates the coordinates of the atoms forming one unit of this pattern.

# The unit cell of the projected pattern is a rectangle with atoms at the four corners
# and one atom at the center. The side lengths of this rectangle are 'a' and 'a * sqrt(2)',
# where 'a' is the lattice constant of the 3D FCC cell.

# Let's set the lattice constant 'a' for our calculation. The exact value does not
# affect the shape, so we choose one for demonstration.
a = 10

# Calculate the dimensions of the rectangular unit cell in the 2D projection.
side_length_1 = a
side_length_2 = a * math.sqrt(2)

print("For an FCC lattice viewed along the [110] direction:")
print("The projected pattern is a centered rectangle.")
print(f"If the lattice constant 'a' is {a}, the rectangle's dimensions are:")
print(f"Side 1 = {side_length_1:.3f}")
print(f"Side 2 = a * sqrt(2) = {side_length_2:.3f}")
print("")

# Calculate the coordinates of the atoms for one rectangular unit.
# We place one corner at the origin (0, 0) for simplicity.
corner_1 = (0.0, 0.0)
corner_2 = (side_length_2, 0.0)
corner_3 = (0.0, side_length_1)
corner_4 = (side_length_2, side_length_1)
center   = (side_length_2 / 2.0, side_length_1 / 2.0)

print("The 2D coordinates for the atoms in one unit cell of the pattern are:")
print(f"Corner 1: ({corner_1[0]:.3f}, {corner_1[1]:.3f})")
print(f"Corner 2: ({corner_2[0]:.3f}, {corner_2[1]:.3f})")
print(f"Corner 3: ({corner_3[0]:.3f}, {corner_3[1]:.3f})")
print(f"Corner 4: ({corner_4[0]:.3f}, {corner_4[1]:.3f})")
print(f"Center  : ({center[0]:.3f}, {center[1]:.3f})")
print("\nBy visually inspecting the options, Image D is the only one that displays this centered rectangular pattern.")