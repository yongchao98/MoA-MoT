import numpy as np

# Step 1 & 2: Define the model and parameters.
# The problem boils down to finding the center of mass of a composite object.
# Through torque analysis on the entire string system, it can be shown that
# for the system to be in equilibrium, the length of the hanging part (L_h)
# must be equal to the radius of the hemisphere (R).
# Let's set the radius R to 1 for a numerical answer.
R = 1.0
L_h = R

# We assume a constant mass per unit length, rho. It will cancel out.
rho = 1.0

# The string consists of two parts:
# Part 1: A quarter-circular arc on the hemisphere.
# Part 2: A vertical line segment hanging down.

# Step 3: Calculate properties of Part 1 (the arc).
# The arc is a quarter circle from theta=0 to theta=pi/2.
# x = R*sin(theta), z = R*cos(theta)
length_arc = R * np.pi / 2
mass_arc = rho * length_arc

# The center of mass of a quarter-circular arc of radius R is at
# (x, z) = (2*R/pi, 2*R/pi) in a coordinate system where the arc is in the first quadrant.
# Our setup matches this, starting from z-axis (theta=0) to x-axis (theta=pi/2).
x_cm_arc = 2 * R / np.pi
z_cm_arc = 2 * R / np.pi

# Step 4: Calculate properties of Part 2 (the hanging segment).
length_hang = L_h
mass_hang = rho * length_hang

# The hanging part starts at the equator, point P, with coordinates (R, 0).
# It hangs vertically downwards for a length L_h = R.
# Its starting point is (R, 0) and its end point is (R, -R).
# The center of mass for this vertical line is at its midpoint.
x_cm_hang = R
z_cm_hang = (0 + (-R)) / 2.0

# Step 5: Calculate the center of mass of the combined system.
total_mass = mass_arc + mass_hang

# The horizontal coordinate (x) of the total center of mass:
x_cm_total = (mass_arc * x_cm_arc + mass_hang * x_cm_hang) / total_mass
x_cm_total_formula = (rho * (R * np.pi / 2) * (2 * R / np.pi) + rho * R * R) / (rho * (R * np.pi / 2) + rho * R)
# Simplifying the formula:
# x_cm_total_formula = (R**2 + R**2) / (R * np.pi / 2 + R)
# x_cm_total_formula = 2 * R**2 / (R * (np.pi/2 + 1))
# x_cm_total_formula = 2 * R / (np.pi/2 + 1)
# x_cm_total_formula = 4 * R / (np.pi + 2)
final_x = 4 * R / (np.pi + 2)

# The vertical coordinate (z) of the total center of mass:
z_cm_total = (mass_arc * z_cm_arc + mass_hang * z_cm_hang) / total_mass
z_cm_total_formula = (rho * (R * np.pi / 2) * (2 * R / np.pi) + rho * R * (-R / 2)) / (rho * (R * np.pi / 2) + rho * R)
# Simplifying the formula:
# z_cm_total_formula = (R**2 - R**2 / 2) / (R * np.pi / 2 + R)
# z_cm_total_formula = (R**2 / 2) / (R * (np.pi/2 + 1))
# z_cm_total_formula = (R / 2) / (np.pi/2 + 1)
# z_cm_total_formula = R / (np.pi + 2)
final_z = R / (np.pi + 2)

# Print the final raw numbers for the coordinates, separated by a comma.
print(f"The equation for the horizontal coordinate is (4 * R) / (pi + 2).")
print(f"For R=1, the horizontal coordinate is 4 / ({np.pi} + 2) = {final_x}")
print(f"The equation for the vertical coordinate is R / (pi + 2).")
print(f"For R=1, the vertical coordinate is 1 / ({np.pi} + 2) = {final_z}")
print(f"\nThe raw numbers for the horizontal and vertical coordinates are:")
print(f"{final_x},{final_z}")
