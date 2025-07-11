import math

# Define the constants based on the problem description
rho = 1000  # Density of water in kg/m^3
g = 10      # Acceleration due to gravity in m/s^2
H = 10      # Depth of the river in meters

# The initial pressure at the bottom is the hydrostatic pressure P = ρ * g * H
initial_pressure = rho * g * H

# The pressure drop due to flow is the dynamic pressure = 1/2 * ρ * v^2
# We need to find the speed 'v' for which the pressure at the bottom becomes 0.
# This occurs when the dynamic pressure equals the initial hydrostatic pressure.

# Set up the equation:
# initial_pressure = 1/2 * ρ * v^2
# ρ * g * H = 1/2 * ρ * v^2

print("We are solving the equation for the flow speed 'v'.")
print("The initial pressure at the bottom (ρ * g * H) must be equal to the dynamic pressure (1/2 * ρ * v^2) for the net pressure to be zero.\n")
print("Equation: ρ * g * H = 1/2 * ρ * v^2")
print(f"Substituting the values: {rho} * {g} * {H} = 1/2 * {rho} * v^2")
print(f"{initial_pressure} = {0.5 * rho} * v^2\n")

print("To solve for v, we can first simplify the equation to v^2 = 2 * g * H")
v_squared = 2 * g * H
print(f"v^2 = 2 * {g} * {H}")
print(f"v^2 = {v_squared}\n")

# Calculate v by taking the square root
v = math.sqrt(v_squared)

print("Now, we find v by taking the square root:")
print(f"v = sqrt({v_squared})")
print(f"The flow speed v is {v:.2f} m/s.")
