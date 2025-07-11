import math

# Based on the physical analysis of the system, the center of mass coordinates
# are given by specific formulas. We assume a radius R=1, as the final
# coordinates are proportional to R.

# Radius of the sphere
R = 1.0
pi = math.pi

# Part 1: The string on the sphere (arc from theta=0 to pi/2)
# Length and effective mass (proportional to length)
L1 = R * pi / 2
M1 = L1
# Center of mass of the arc
X1 = 2 * R / pi
Z1 = 2 * R / pi

# Part 2: The string hanging vertically
# Length and effective mass
L2 = R
M2 = L2
# Center of mass of the hanging segment
X2 = R
Z2 = -R / 2

# Total mass
M_total = M1 + M2

# Calculate the horizontal (X) and vertical (Z) coordinates of the total center of mass
# X_cm = (M1*X1 + M2*X2) / M_total = ( (R*pi/2)*(2*R/pi) + R*R ) / (R*pi/2 + R)
#      = (R^2 + R^2) / (R*(pi/2 + 1)) = 2*R^2 / (R*(pi+2)/2) = 4*R / (pi+2)
X_cm = 4 * R / (pi + 2)

# Z_cm = (M1*Z1 + M2*Z2) / M_total = ( (R*pi/2)*(2*R/pi) + R*(-R/2) ) / (R*pi/2 + R)
#      = (R^2 - R^2/2) / (R*(pi/2 + 1)) = (R^2/2) / (R*(pi+2)/2) = R / (pi+2)
Z_cm = R / (pi + 2)

# The question asks for the raw number of the horizontal and vertical coordinates
# separated by a comma. We print these calculated values.
print(f"{X_cm},{Z_cm}")