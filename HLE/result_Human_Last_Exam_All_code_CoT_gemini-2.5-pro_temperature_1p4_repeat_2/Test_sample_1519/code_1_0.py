import math

# Define the given constants for the problem
g = 10.0  # Acceleration due to gravity in m/s^2
H = 10.0  # Depth of the river in meters
# The density of water (ρ = 1000 kg/m^3) is not needed as it cancels out.

print("According to the problem, the pressure at the bottom of the river becomes zero.")
print("This occurs when the initial hydrostatic pressure (ρ*g*H) is completely converted into dynamic pressure (1/2*ρ*v^2).")
print("The equation is: ρ * g * H = 1/2 * ρ * v^2")
print("By simplifying, we can solve for the velocity v: v = sqrt(2 * g * H)")
print("\nNow, we substitute the given values into the equation:")

# Calculate v^2
v_squared = 2 * g * H
# Calculate v
v = math.sqrt(v_squared)

# Output the final equation with the numbers and the result.
print(f"v = sqrt(2 * {g} * {H})")
print(f"v = sqrt({v_squared})")
print(f"v = {v} m/s")