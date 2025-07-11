import math

# Step 1 & 2: Define the geometry and parameters.
# We can set R=1 and rho=1 as they will cancel out when calculating the
# dimensionless coordinates of the center of mass.
R = 1.0

# The path of the string is a quarter-circle arc on the sphere (from theta=0 to pi/2)
# followed by a vertical segment hanging down.

# Step 3: The length of the hanging part 'h' is equal to the radius 'R'.
h = R

# Step 4: Define the properties of the two parts of the string.
# Part 1: The arc
m_arc = (math.pi * R) / 2
# The center of mass of a quarter-circle arc from theta=0 to pi/2 in the x-z plane
# is at (2R/pi, 2R/pi).
x_cm_arc = 2 * R / math.pi
z_cm_arc = 2 * R / math.pi

# Part 2: The hanging vertical segment
m_hang = h
# The segment hangs from (R, 0, 0) down to (R, 0, -h). Its center of mass is at its midpoint.
x_cm_hang = R
z_cm_hang = -h / 2

# Step 5: Calculate the total mass and the coordinates of the overall center of mass.
total_mass = m_arc + m_hang

# The horizontal coordinate of the center of mass (X_cm)
# We use the formula: X_cm = (m_arc * x_cm_arc + m_hang * x_cm_hang) / total_mass
x_cm = (m_arc * x_cm_arc + m_hang * x_cm_hang) / total_mass

# The vertical coordinate of the center of mass (Z_cm)
# We use the formula: Z_cm = (m_arc * z_cm_arc + m_hang * z_cm_hang) / total_mass
z_cm = (m_arc * z_cm_arc + m_hang * z_cm_hang) / total_mass

# Analytically, the coordinates are:
# X_cm/R = 4 / (2 + pi)
# Z_cm/R = 1 / (2 + pi)

# Step 6: Print the raw numerical values of the coordinates.
# The horizontal coordinate is x_cm and the vertical coordinate is z_cm.
# The question "without considering z-axis coordinates" is interpreted as a request
# for the (horizontal, vertical) pair in the plane of the string, ignoring the third dimension (y-axis).
horizontal_coordinate = 4 / (2 + math.pi)
vertical_coordinate = 1 / (2 + math.pi)

print(f"{horizontal_coordinate},{vertical_coordinate}")