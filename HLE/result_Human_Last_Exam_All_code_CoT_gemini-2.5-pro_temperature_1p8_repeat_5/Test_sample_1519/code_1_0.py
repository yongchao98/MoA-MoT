import math

# Given parameters
rho = 1000  # density of water in kg/m^3
g = 10      # acceleration due to gravity in m/s^2
H = 10      # depth of the river in meters

# The problem is to find the speed v at which the pressure at the bottom becomes 0.
# The initial hydrostatic pressure (P_rest) is counteracted by the dynamic pressure (P_dynamic).
# P_new = P_rest - P_dynamic
# 0 = (rho * g * H) - (1/2 * rho * v^2)

# Calculate the speed v
# v = sqrt(2 * g * H)
v = math.sqrt(2 * g * H)

# Print the final equation with each number, as requested.
print("The final equation for the pressure at the bottom to be zero is:")
print(f"0 = ({rho} * {g} * {H}) - (1/2 * {rho} * v^2)")

# Print the solution
print("\nSolving for v:")
print(f"v = sqrt(2 * {g} * {H})")
print(f"v = {v} m/s")

# The final answer
print(f'<<<{v}>>>')