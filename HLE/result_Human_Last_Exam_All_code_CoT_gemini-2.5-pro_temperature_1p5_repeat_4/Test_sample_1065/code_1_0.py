import math

# This script calculates the coordinates of the center of mass of the candy string.
# Based on the physical analysis, the string forms a circular arc on the pumpkin's surface,
# spanning from an angle of 0 to pi/3 from the vertical axis.
#
# The radius of the pumpkin is R. We are calculating the numerical coefficients
# for the center of mass coordinates (X_cm, Z_cm).
#
# The derived symbolic formulas for the coefficients are:
# Horizontal coefficient: 3 / (2 * pi)
# Vertical coefficient: (3 * sqrt(3)) / (2 * pi)
#
# The numbers in the final equation for the horizontal coordinate are 3, 2, and pi.
print("Calculating the horizontal coordinate coefficient from the expression: 3 / (2 * pi)")
horizontal_coeff = 3 / (2 * math.pi)

# The numbers in the final equation for the vertical coordinate are 3, sqrt(3), 2, and pi.
print("Calculating the vertical coordinate coefficient from the expression: (3 * sqrt(3)) / (2 * pi)")
vertical_coeff = (3 * math.sqrt(3)) / (2 * math.pi)

# The problem asks for the raw numbers of the horizontal and vertical coordinates, separated by a comma.
print("\nThe final numerical coefficients are:")
print(f"{horizontal_coeff},{vertical_coeff}")
