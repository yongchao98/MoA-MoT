import math

# The problem asks for the horizontal and vertical coordinates of the center of mass.
# Based on our physical model, the coordinates are given by the formulas:
# X_cm = 4*R / (pi + 2)
# Z_cm = R / (pi + 2)
# We will calculate the numerical coefficients by assuming R=1.

# Calculate the coefficient for the horizontal coordinate (x)
x_coefficient = 4 / (math.pi + 2)

# Calculate the coefficient for the vertical coordinate (z)
z_coefficient = 1 / (math.pi + 2)

# The problem asks for the raw numbers of the horizontal and vertical coordinates separated by a comma.
# The following code prints these two values.
print(f"{x_coefficient},{z_coefficient}")
