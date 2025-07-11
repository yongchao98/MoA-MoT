import math

# This script calculates the horizontal and vertical coordinates of the center of mass
# for the string of candies. The coordinates are expressed as coefficients of the radius R.

# The derivation leads to the following formulas for the center of mass (X_cm, Z_cm):
# X_cm = (4 * R) / (pi + 2)
# Z_cm = (R) / (pi + 2)

# The problem asks for the raw numerical coefficients for these coordinates.
pi = math.pi

# --- Final Equation Numbers ---
# To satisfy the instruction "output each number in the final equation", we show
# the derivation of the final expressions. The CM calculation combines the arc
# segment and the vertical segment. Let R=1 and mass density rho=1 for simplicity
# as these terms cancel out in the final coefficient calculation.

# Part 1: Quarter-circle arc
m_arc = 1 * (pi/2)  # Mass = rho * length
x_arc = 2 / pi      # x-coordinate of CM of arc
z_arc = 2 / pi      # z-coordinate of CM of arc

# Part 2: Vertical hanging segment
m_vert = 1          # Mass = rho * length (where length = R=1)
x_vert = 1          # x-coordinate of CM of vertical part
z_vert = -1 / 2     # z-coordinate of CM of vertical part

# The "final equation" combines these numbers:
# X_cm = (m_arc*x_arc + m_vert*x_vert) / (m_arc + m_vert)
# Z_cm = (m_arc*z_arc + m_vert*z_vert) / (m_arc + m_vert)
print("The numbers entering the final calculation (assuming R=1, rho=1) are:")
print(f"Arc mass: {m_arc}")
print(f"Arc CM x-coordinate: {x_arc}")
print(f"Arc CM z-coordinate: {z_arc}")
print(f"Vertical mass: {m_vert}")
print(f"Vertical CM x-coordinate: {x_vert}")
print(f"Vertical CM z-coordinate: {z_vert}")
print("-" * 20)

# Calculate the final coefficients
horizontal_coeff = 4 / (pi + 2)
vertical_coeff = 1 / (pi + 2)

# Print the final result as comma-separated raw numbers
print("The horizontal and vertical coordinates' raw numerical coefficients are:")
print(f"{horizontal_coeff},{vertical_coeff}")