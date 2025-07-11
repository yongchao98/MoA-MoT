import math

# Define the constants based on the problem description
rho = 1000  # Density of water in kg/m^3
g = 10      # Acceleration due to gravity in m/s^2
H = 10      # Depth of the river in meters

# The condition is that the static pressure equals the dynamic pressure
# for the final pressure to be zero.
# rho * g * H = 0.5 * rho * v^2
# We can simplify this to g * H = 0.5 * v^2
# Solving for v: v = sqrt(2 * g * H)

# Calculate the flow speed v
v_squared = 2 * g * H
v = math.sqrt(v_squared)

# Print the equation with the values substituted
print("The relationship is given by: rho * g * H = (1/2) * rho * v^2")
print("This simplifies to: g * H = 0.5 * v^2")
print("Solving for v, we get: v = sqrt(2 * g * H)")
print(f"v = sqrt(2 * {g} * {H})")
print(f"v = sqrt({v_squared})")
print(f"The flow speed v is: {v} m/s")

# Final answer in the specified format
print(f"\n<<<{v}>>>")