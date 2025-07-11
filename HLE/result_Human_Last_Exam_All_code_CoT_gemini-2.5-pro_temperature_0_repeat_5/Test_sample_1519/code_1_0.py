import math

# Define the given physical constants and parameters
rho = 1000  # Density of water in kg/m^3
g = 10      # Acceleration due to gravity in m/s^2
H = 10      # Depth of the river in meters

# According to Bernoulli's principle, the initial hydrostatic pressure
# at the bottom of the river is converted into dynamic pressure when the water flows.
# We want to find the speed 'v' where the pressure at the bottom becomes zero.
# The equation is: P_initial = P_dynamic
# rho * g * H = (1/2) * rho * v^2
#
# We can simplify this by canceling rho from both sides and solving for v:
# g * H = (1/2) * v^2
# v^2 = 2 * g * H
# v = sqrt(2 * g * H)

# Calculate the value of v
v = math.sqrt(2 * g * H)

# Print the final equation with the substituted values and the result
print("To find the flow speed v, we use the equation: v = sqrt(2 * g * H)")
print("Substituting the values:")
print(f"v = sqrt(2 * {g} * {H})")
print(f"The required flow speed is {v} m/s.")