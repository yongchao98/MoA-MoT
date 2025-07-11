import math

# Define the given constants
rho = 1000  # Density of water in kg/m^3
g = 10      # Acceleration due to gravity in m/s^2
H = 10      # Depth of the river in meters

# The equation for the gauge pressure at the bottom is:
# P_final = (rho * g * H) - (0.5 * rho * v^2)

# We set P_final to 0 and solve for v.
# 0 = (rho * g * H) - (0.5 * rho * v^2)
# v = sqrt(2 * g * H)

# Calculate the speed v
v_squared = 2 * g * H
v = math.sqrt(v_squared)

# Print the final equation with the numbers plugged in
print("The final equation to solve for the speed 'v' where the pressure becomes 0 is:")
print(f"0 = ({rho} * {g} * {H}) - (0.5 * {rho} * v^2)")
print(f"This simplifies to v = sqrt(2 * {g} * {H})")
print(f"v = sqrt({v_squared})")
print(f"\nThe flow speed v at which the pressure at the bottom decreases to zero is: {v} m/s")