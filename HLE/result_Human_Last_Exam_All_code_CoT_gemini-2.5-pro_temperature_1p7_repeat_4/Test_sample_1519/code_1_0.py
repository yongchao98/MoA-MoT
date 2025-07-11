import math

# Define the given physical constants
rho = 1000  # Density of water in kg/m^3
g = 10      # Acceleration due to gravity in m/s^2
H = 10      # Depth of the river in meters

# The initial pressure at the bottom (P_static) is given by rho * g * H.
# When the water flows, the pressure at the bottom (P_new) is reduced by the dynamic pressure (0.5 * rho * v^2).
# P_new = (rho * g * H) - (0.5 * rho * v^2)

# We need to find the speed 'v' for which the pressure at the bottom becomes 0.
# So, we solve the equation: 0 = (rho * g * H) - (0.5 * rho * v^2)

print("We are solving for the flow speed 'v' where the pressure at the bottom of the river becomes 0.")
print("The equation to solve is:")
print(f"0 = ({rho} * {g} * {H}) - (0.5 * {rho} * v^2)\n")

# To solve for v, we can rearrange the equation:
# 0.5 * rho * v^2 = rho * g * H
# v^2 = 2 * g * H
v_squared = 2 * g * H
v = math.sqrt(v_squared)

# Print the calculation steps
p_static = rho * g * H
print(f"First, calculate the static pressure component: {rho} * {g} * {H} = {p_static} N/m^2")
print(f"Our equation becomes: 0 = {p_static} - (0.5 * {rho} * v^2)")
print(f"Rearranging for v^2: v^2 = {p_static} / (0.5 * {rho})")
print(f"v^2 = {v_squared}")
print(f"Solving for v: v = sqrt({v_squared})")

print(f"\nThe required flow speed v is approximately {v:.2f} m/s.")