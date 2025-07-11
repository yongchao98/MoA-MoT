import math

# Step 1: Define the constants.
# We are looking for the coordinates as coefficients of the radius R.
# We can set R=1 for this calculation.
# pi is the mathematical constant.
R = 1
pi = math.pi

# Step 2: Define the properties of the two parts of the string.
# Part 1: The arc on the hemisphere (a quarter circle).
# The length of the arc is L_arc = R * (pi / 2).
# The mass of the arc is m_arc = rho * L_arc = rho * R * pi / 2.
# The center of mass of a quarter-circle arc of radius R is at (2*R/pi, 2*R/pi).
# Let's set rho=1 as it will cancel out.
m_arc = R * pi / 2
x_arc = 2 * R / pi
z_arc = 2 * R / pi

# Part 2: The hanging vertical line.
# The length of the hanging part is L_line = R.
# The mass of the line is m_line = rho * L_line = rho * R.
# The line hangs from (R, 0) down to (R, -R).
# Its center of mass is at (R, -R/2).
m_line = R
x_line = R
z_line = -R / 2

# Step 3: Calculate the total mass.
# The total mass M = m_arc + m_line.
M = m_arc + m_line

# Step 4: Calculate the coordinates of the overall Center of Mass (CM).
# The formula for the CM of a composite body is:
# X_cm = (m_arc * x_arc + m_line * x_line) / M
# Z_cm = (m_arc * z_arc + m_line * z_line) / M

# Horizontal coordinate (X_cm)
X_cm_numerator = m_arc * x_arc + m_line * x_line
X_cm = X_cm_numerator / M

# Vertical coordinate (Z_cm)
Z_cm_numerator = m_arc * z_arc + m_line * z_line
Z_cm = Z_cm_numerator / M

# Step 5: Print the results and the formulas.
# The analytical expressions for the coordinates are:
# X_cm = ( (rho*R*pi/2)*(2*R/pi) + (rho*R)*R ) / (rho*R*pi/2 + rho*R)
#      = ( R^2 + R^2 ) / (R*pi/2 + R) = 2*R^2 / (R*(pi/2 + 1))
#      = 2*R / (pi/2 + 1) = 4*R / (pi + 2)
#
# Z_cm = ( (rho*R*pi/2)*(2*R/pi) + (rho*R)*(-R/2) ) / (rho*R*pi/2 + rho*R)
#      = ( R^2 - R^2/2 ) / (R*pi/2 + R) = (R^2/2) / (R*(pi/2 + 1))
#      = (R/2) / (pi/2 + 1) = R / (pi + 2)

print("The final equations for the center of mass coordinates (X_cm, Z_cm) are:")
print("Horizontal coordinate X_cm = (4 * R) / (pi + 2)")
print("Vertical coordinate Z_cm = (1 * R) / (pi + 2)")
print("\nThe raw numerical coefficients (for R=1) are:")

# The final output requires printing the raw numbers separated by a comma.
print(f"{X_cm},{Z_cm}")