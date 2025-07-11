import math

# Given values
rho = 1000  # Density of water in kg/m^3
g = 10      # Acceleration due to gravity in m/s^2
H = 10      # Depth of the river in meters

# The pressure at the bottom becomes zero when the initial hydrostatic pressure
# is completely balanced by the dynamic pressure from the flow.
#
# Hydrostatic Pressure = Dynamic Pressure
# ρ * g * H = (1/2) * ρ * v^2

print("To find the flow speed v where the pressure at the bottom is zero, we set the hydrostatic pressure equal to the dynamic pressure.")
print("Equation: ρ * g * H = 1/2 * ρ * v^2\n")

# Print the equation with the given values
print(f"Substituting the values:")
print(f"{rho} * {g} * {H} = 0.5 * {rho} * v^2\n")

# Calculate the left side
hydrostatic_pressure = rho * g * H
print(f"Calculating the left side:")
print(f"{int(hydrostatic_pressure)} = 0.5 * {rho} * v^2\n")

# The equation can be simplified by dividing both sides by ρ
print("Simplifying the equation by removing ρ from both sides:")
print(f"{g} * {H} = 0.5 * v^2")
print(f"{g * H} = 0.5 * v^2\n")

# Solve for v^2
v_squared = 2 * g * H
print("Solving for v^2:")
print(f"v^2 = 2 * {g * H}")
print(f"v^2 = {v_squared}\n")

# Solve for v
v = math.sqrt(v_squared)
print("Solving for v by taking the square root:")
print(f"v = sqrt({v_squared})")
print(f"v = {v}\n")

print(f"Therefore, the flow speed required for the pressure at the bottom to decrease to zero is {v:.2f} m/s.")