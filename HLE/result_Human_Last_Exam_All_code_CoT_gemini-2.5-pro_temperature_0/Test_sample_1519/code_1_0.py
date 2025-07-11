import math

# Define the given constants
g = 10  # Acceleration due to gravity in m/s^2
H = 10  # Depth of the river in meters

# The equation to solve for v is v = sqrt(2 * g * H)
# We first calculate the value inside the square root
v_squared = 2 * g * H

# Then we take the square root to find v
v = math.sqrt(v_squared)

# Print the steps of the final calculation
print("To find the flow speed 'v' where the pressure at the bottom is zero, we solve the equation:")
print("v = sqrt(2 * g * H)")
print("Substituting the values:")
print(f"v = sqrt(2 * {g} * {H})")
print(f"v = sqrt({v_squared})")
print(f"The required flow speed is: {v} m/s")