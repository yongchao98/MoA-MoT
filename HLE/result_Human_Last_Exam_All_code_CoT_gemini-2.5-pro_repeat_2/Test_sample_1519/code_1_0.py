import math

# Step 1: Define the given parameters
rho = 1000  # Density of water in kg/m^3
g = 10      # Acceleration due to gravity in m/s^2
H = 10      # Depth of the river in meters

# Step 2: Calculate the initial pressure at the bottom when the water is at rest
P_rest = rho * g * H

# Step 3: Set up the equation for the final state
# The pressure in the flowing water is P_flow = P_rest - (1/2) * rho * v^2
# We need to find the speed v when P_flow = 0.
# The equation becomes: 0 = P_rest - (1/2) * rho * v^2

# Step 4: Rearrange the equation to solve for v
# (1/2) * rho * v^2 = P_rest
# (1/2) * rho * v^2 = rho * g * H
# v^2 = 2 * g * H
v_squared = 2 * g * H
v = math.sqrt(v_squared)

# Step 5: Print the steps with the numerical values
print("The goal is to find the flow speed v when the pressure at the bottom becomes 0.")
print("The initial pressure at rest is P_rest = ρ * g * H.")
print(f"P_rest = {rho} * {g} * {H} = {P_rest} N/m^2\n")

print("The pressure when flowing, P_flow, is set to 0:")
print("P_flow = P_rest - (1/2) * ρ * v^2")
print(f"0 = {P_rest} - (1/2) * {rho} * v^2\n")

print("Solving for v:")
print(f"(1/2) * {rho} * v^2 = {P_rest}")
print(f"{0.5 * rho} * v^2 = {P_rest}")
print(f"v^2 = {P_rest} / {0.5 * rho}")
print(f"v^2 = {v_squared}")
print(f"v = sqrt({v_squared})\n")

print(f"The required flow speed is {v:.2f} m/s.")

print(f"\n<<<{v:.2f}>>>")