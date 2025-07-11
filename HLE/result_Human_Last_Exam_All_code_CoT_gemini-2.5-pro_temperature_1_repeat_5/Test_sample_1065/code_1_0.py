import math

# Plan:
# 1. Model the system as a string on a hemisphere of radius R.
# 2. The string is in equilibrium, with the part on the surface forming a quarter-circle arc.
# 3. Assume the mass of the hanging part (M_hang) is 1/4 the mass of the arc part (M_arc).
# 4. Calculate the center of mass (CM) for the arc and the hanging part separately.
# 5. Combine the CMs to find the total CM of the string.
# 6. Provide the final coordinates as raw numbers (in units of R).

# For simplicity, we can assume R=1 and the mass density rho=1.
# The results will be the coordinates of the CM divided by R.
R = 1.0

# --- Part 1: The arc on the hemisphere ---
# The arc is a quarter circle, so its length is (2 * pi * R) / 4.
length_arc = R * math.pi / 2.0
mass_arc = length_arc  # Since rho = 1

# The CM of a quarter-circle arc in the first quadrant is at (2R/pi, 2R/pi).
x_arc = 2 * R / math.pi
z_arc = 2 * R / math.pi

# --- Part 2: The hanging string ---
# Using the assumption from the plan.
mass_hang = 0.25 * mass_arc
length_hang = mass_hang # Since rho = 1

# The hanging part starts at (R, 0) and hangs vertically down.
# Its CM is at its midpoint.
x_hang = R
z_hang = -length_hang / 2.0

# --- Part 3: Combined Center of Mass ---
# Total mass of the string.
mass_total = mass_arc + mass_hang

# Calculate the weighted average for the x and z coordinates of the CM.
# x_cm = (x_arc * mass_arc + x_hang * mass_hang) / mass_total
# z_cm = (z_arc * mass_arc + z_hang * mass_hang) / mass_total
x_cm = (x_arc * mass_arc + x_hang * mass_hang) / mass_total
z_cm = (z_arc * mass_arc + z_hang * mass_hang) / mass_total

# The problem asks for the raw number of the horizontal and vertical coordinates.
# The calculation gives these coordinates in units of R.
print(f"{x_cm},{z_cm}")
