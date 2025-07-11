import math

# Define the given physical constants
g = 10  # Acceleration due to gravity in m/s^2
H = 10  # Depth of the river in meters

# According to Bernoulli's principle, for the pressure at the bottom to be zero,
# the initial hydrostatic pressure (rho * g * H) must equal the dynamic pressure (1/2 * rho * v^2).
# rho * g * H = 1/2 * rho * v^2
# g * H = 1/2 * v^2
# v^2 = 2 * g * H
# v = sqrt(2 * g * H)

# Calculate the value of v^2
v_squared = 2 * g * H

# Calculate the final speed v
v = math.sqrt(v_squared)

# Print the calculation steps
print("To find the flow speed v, we set the initial hydrostatic pressure equal to the dynamic pressure:")
print("rho * g * H = 1/2 * rho * v^2")
print("This simplifies to v = sqrt(2 * g * H)")
print("\nPlugging in the values:")
print(f"v = sqrt(2 * {g} * {H})")
print(f"v = sqrt({v_squared})")
print(f"\nThe flow speed at which the pressure becomes zero is: {v:.2f} m/s")