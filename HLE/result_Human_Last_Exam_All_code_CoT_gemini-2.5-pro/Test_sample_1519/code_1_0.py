import math

# Define the given physical constants
rho = 1000  # Density of water in kg/m^3
g = 10      # Acceleration due to gravity in m/s^2
H = 10      # Depth of the river in meters

# The problem is to find the speed 'v' where the gauge pressure at the bottom becomes zero.
# This occurs when the static pressure from the water column is balanced by the dynamic pressure from the flow.

# The static pressure at the bottom is P_static = ρ * g * H
# The dynamic pressure is P_dynamic = (1/2) * ρ * v^2

# We set them equal to each other: P_static = P_dynamic
print("To find the flow speed 'v' where the pressure at the bottom becomes zero, we set the static pressure equal to the dynamic pressure:")
print("Static Pressure = Dynamic Pressure")
print("ρ * g * H = (1/2) * ρ * v^2")
print("\nWe can solve for 'v'. First, let's substitute the values into the equation:")

# We need to print the equation with all the numbers.
# rho * g * H = (1/2) * rho * v^2
# 1000 * 10 * 10 = (1/2) * 1000 * v^2
print(f"{rho} * {g} * {H} = (1/2) * {rho} * v^2")

# Calculate the left side
static_pressure_value = rho * g * H
print(f"{static_pressure_value} = {0.5 * rho} * v^2")

# Isolate v^2
v_squared = (rho * g * H) / (0.5 * rho)
print(f"v^2 = {static_pressure_value} / {0.5 * rho}")
print(f"v^2 = {v_squared}")


# Calculate v
v = math.sqrt(v_squared)
print(f"v = sqrt({v_squared})")
print(f"\nThe required flow speed is v = {v} m/s.")

print("<<<14.142135623730951>>>")