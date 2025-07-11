import math

# Define the constants from the problem description
rho = 1000  # Density of water in kg/m^3
g = 10      # Acceleration due to gravity in m/s^2
H = 10      # Depth of the river in meters

# Calculate the initial gauge pressure at the bottom when the water is at rest
initial_pressure = rho * g * H

# To find the speed 'v' at which the pressure at the bottom becomes zero,
# we set the initial static pressure equal to the dynamic pressure.
# The equation is: P_initial = (1/2) * rho * v^2
# This simplifies to: rho * g * H = (1/2) * rho * v^2
# Which further simplifies to: v = sqrt(2 * g * H)

# Solve for the velocity v
v = math.sqrt(2 * g * H)

# Print the step-by-step equation and the final answer
print("The initial gauge pressure at the bottom of the river is calculated as:")
print(f"P = ρ * g * H")
print(f"P = {rho} * {g} * {H} = {int(initial_pressure)} N/m^2")
print("\nTo find the speed 'v' where this pressure decreases to zero, we set the initial pressure equal to the dynamic pressure term:")
print(f"ρ * g * H = (1/2) * ρ * v^2")
print(f"{int(initial_pressure)} = (1/2) * {rho} * v^2")
print("\nSolving for v:")
print(f"v = sqrt(2 * g * H)")
print(f"v = sqrt(2 * {g} * {H})")
print(f"v = sqrt({2 * g * H})")
print(f"\nThe required flow speed is v = {v:.2f} m/s.")
<<<14.14>>>