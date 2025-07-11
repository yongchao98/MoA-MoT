import math

# Define the constants based on the problem description
g = 10  # acceleration due to gravity in m/s^2
H = 10  # depth of the river in meters

# The equation to find the speed v where the pressure at the bottom becomes zero is
# derived from Bernoulli's principle:
# Initial Pressure = Pressure drop due to flow
# rho * g * H = 1/2 * rho * v^2
# This simplifies to v = sqrt(2 * g * H)

# Calculate the speed v
v = math.sqrt(2 * g * H)

# Print the final equation with the numerical values and the result.
# The user wants to see each number in the final equation.
print(f"To find the flow speed 'v' where the pressure becomes zero, we solve the equation:")
print(f"v = sqrt(2 * g * H)")
print(f"v = sqrt(2 * {g} * {H})")
print(f"v = {v:.2f} m/s")