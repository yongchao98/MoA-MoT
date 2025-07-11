import math

# Step 1: Define the given physical quantities.
rho = 1000  # Density of water in kg/m^3
g = 10      # Acceleration due to gravity in m/s^2
H = 10      # Depth of the river in meters

# Step 2: Calculate the initial static pressure at the bottom.
P_static = rho * g * H

# Step 3 & 4: The problem states that the final pressure is zero. 
# This means the initial static pressure is equal to the dynamic pressure (1/2 * rho * v^2).
# We need to find the velocity 'v' that satisfies this condition.

# Step 5 & 6: Set up and solve the equation for v.
# P_static = (1/2) * rho * v^2
# Rearranging for v^2: v^2 = (2 * P_static) / rho
# Since P_static = rho * g * H, we can substitute it in:
# v^2 = (2 * rho * g * H) / rho
# v^2 = 2 * g * H
v = math.sqrt(2 * g * H)

# Print the final equation with the numerical values as requested.
print(f"The equation to solve is based on the initial pressure at the bottom being converted to dynamic pressure:")
print(f"{int(P_static)} = 1/2 * {int(rho)} * v^2")
print(f"\nSolving for the flow speed v, we get:")
print(f"v = sqrt(2 * {g} * {H}) = {v} m/s")
