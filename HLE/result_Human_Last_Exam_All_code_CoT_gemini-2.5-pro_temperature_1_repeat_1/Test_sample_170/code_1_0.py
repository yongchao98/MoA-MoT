import math

# Define the anisotropic ratio
AR = 0.1

# --- Part 1: Calculate the optimal orientation angle (theta) ---

# The optimal orientation is found by solving the equation:
# cos(2 * theta) = (1 - AR) / (1 + AR)
cos_2_theta = (1 - AR) / (1 + AR)

# Solve for 2*theta in radians
two_theta_rad = math.acos(cos_2_theta)

# Convert to degrees and solve for theta
theta_deg = math.degrees(two_theta_rad) / 2

# --- Part 2: Calculate the smallest pressure gradient angle (phi_y) ---

# At the optimal orientation, the tangent of the pressure gradient angle is given by:
# tan(phi_y) = -2 * sqrt(AR) / (1 - AR)
tan_phi_y = -2 * math.sqrt(AR) / (1 - AR)

# Solve for the angle phi_y in radians
phi_y_rad = math.atan(tan_phi_y)

# Convert to degrees and take the absolute value for the final angle
abs_phi_y_deg = abs(math.degrees(phi_y_rad))


# --- Print the results ---

print("To achieve the smallest angle for the pressure gradient, the textile must be oriented at a specific angle.")
print("The orientation angle (theta) of the textile's major permeability axis relative to the flow direction is found by solving:")
print(f"cos(2 * theta) = (1 - {AR}) / (1 + {AR}) = {cos_2_theta:.4f}")
print(f"This gives an orientation angle theta = {theta_deg:.2f} degrees.")

print("\nAt this orientation, the pressure gradient forms an angle (phi_y) with the direction perpendicular to the flow.")
print("This angle is the smallest possible and is calculated from:")
print(f"tan(phi_y) = -2 * sqrt({AR}) / (1 - {AR}) = {tan_phi_y:.4f}")
print(f"The smallest angle is the magnitude of phi_y, which is |phi_y| = {abs_phi_y_deg:.2f} degrees.")
