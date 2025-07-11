import math

# Define the constants
rho = 1000  # Density of water in kg/m^3
g = 10      # Acceleration due to gravity in m/s^2
H = 10      # Depth of the river in meters

# The initial pressure at the bottom when the river is at rest.
P_static = rho * g * H

# The equation to solve is derived from Bernoulli's principle:
# P_flowing + 1/2 * rho * v^2 = P_static
# We want to find v when P_flowing = 0, so the equation becomes:
# 1/2 * rho * v^2 = rho * g * H

# We can solve for v: v = sqrt(2 * g * H)
v_squared = 2 * g * H
v = math.sqrt(v_squared)

print("According to Bernoulli's principle, the static pressure at the bottom is converted into dynamic pressure as the water flows.")
print("We need to find the speed 'v' where the dynamic pressure equals the initial static pressure.")
print("The equation is: 1/2 * rho * v^2 = rho * g * H")
print(f"Plugging in the values: 1/2 * {rho} * v^2 = {rho} * {g} * {H}")
print(f"This simplifies to: {0.5 * rho} * v^2 = {P_static}")
print(f"Solving for v^2: v^2 = {P_static} / {0.5 * rho} = {v_squared}")
print(f"Solving for v: v = sqrt({v_squared})")
print(f"\nThe required flow speed v is: {v:.2f} m/s")

# The final equation with the calculated speed v is:
# 1/2 * 1000 kg/m^3 * (14.14 m/s)^2 = 1000 kg/m^3 * 10 m/s^2 * 10 m
# 100000 N/m^2 = 100000 N/m^2
print("\nFinal equation check with the calculated speed:")
print(f"0 + 1/2 * {rho} * ({v:.2f})^2 = {0.5 * rho * v**2:.0f}")
print(f"And the initial static pressure was: {rho} * {g} * {H} = {P_static}")