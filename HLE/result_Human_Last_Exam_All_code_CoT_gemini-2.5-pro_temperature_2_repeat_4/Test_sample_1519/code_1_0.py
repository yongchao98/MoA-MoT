import math

# Define the constants based on the problem statement
rho = 1000  # Density of water in kg/m^3
g = 10      # Acceleration due to gravity in m/s^2
H = 10      # Depth of the river in meters

# The target final pressure at the bottom
P_final = 0

# The initial pressure is P_initial = rho * g * H
P_initial = rho * g * H

# The problem states that the final pressure is reduced from the initial pressure
# by the dynamic pressure component (1/2 * rho * v^2).
# The equation to solve is: P_final = (rho * g * H) - (1/2 * rho * v^2)

print("The equation relating the initial hydrostatic pressure and the final pressure with the flow speed 'v' is:")
print(f"Final Pressure = (ρ * g * H) - (1/2 * ρ * v^2)")
print("We want to find the speed 'v' for which the Final Pressure is 0. Here is the equation with the given values:")
print(f"{P_final} = ({rho} * {g} * {H}) - (1/2 * {rho} * v^2)\n")

# To solve for v, we can first simplify the equation:
# 1/2 * rho * v^2 = rho * g * H
# v^2 = 2 * g * H
# v = sqrt(2 * g * H)

# Perform the calculation
v_squared = 2 * g * H
v = math.sqrt(v_squared)

print(f"Solving for v:")
print(f"v = sqrt(2 * {g} * {H})")
print(f"v = sqrt({v_squared})")
print(f"\nThe required flow speed 'v' for the pressure at the bottom to be zero is {v:.2f} m/s.")