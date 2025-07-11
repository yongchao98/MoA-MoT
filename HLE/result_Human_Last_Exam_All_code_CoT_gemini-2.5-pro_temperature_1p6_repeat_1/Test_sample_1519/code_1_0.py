import math

# Define the given values
g = 10  # acceleration due to gravity in m/s^2
H = 10  # depth of the river in meters

# The problem is to find the speed v where the pressure at the bottom becomes zero.
# This occurs when the initial hydrostatic pressure (rho * g * H) is completely
# counteracted by the dynamic pressure term (1/2 * rho * v^2).
#
# Equation:
# rho * g * H = (1/2) * rho * v^2
#
# We can simplify by dividing by rho:
# g * H = (1/2) * v^2
#
# Solving for v:
# v^2 = 2 * g * H
# v = sqrt(2 * g * H)

# Calculate the speed v
v_squared = 2 * g * H
v = math.sqrt(v_squared)

print("The equation to solve for the flow speed v is derived by equating the hydrostatic pressure to the dynamic pressure:")
print("v = sqrt(2 * g * H)")
print("\nSubstituting the given values:")
print(f"v = sqrt(2 * {g} * {H})")
print(f"v = sqrt({v_squared})")
print(f"The flow speed v is {v} m/s.")
