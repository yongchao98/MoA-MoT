import math

# Define the given constants
rho = 1000  # Density of water in kg/m^3
g = 10      # Acceleration due to gravity in m/s^2
H = 10      # Depth of the river in meters

# The equation is rho * g * H = 0.5 * rho * v^2
# We solve for v: v = sqrt(2 * g * H)

# Calculate the square of the velocity
v_squared = 2 * g * H

# Calculate the velocity
v = math.sqrt(v_squared)

# The final equation with substituted numbers
# We need to calculate each side of the equation rho*g*H = 1/2*rho*v^2
left_side = rho * g * H
right_side_coeff = 0.5 * rho

print("The pressure at the bottom becomes zero when the initial static pressure equals the dynamic pressure.")
print("The equation is: rho * g * H = (1/2) * rho * v^2")
print("Substituting the values:")
print(f"{rho} * {g} * {H} = 0.5 * {rho} * v^2")
print(f"{left_side} = {right_side_coeff} * v^2")
print(f"v^2 = {left_side} / {right_side_coeff}")
print(f"v^2 = {v_squared}")
print(f"v = sqrt({v_squared})")
print(f"The required flow speed is: {v:.2f} m/s")

# The final answer in the requested format
# <<<14.14>>>