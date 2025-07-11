import math

# Define the physical constants given in the problem
rho = 1000  # Density of water in kg/m^3
g = 10      # Acceleration due to gravity in m/s^2
H = 10      # Depth of the river in meters

# The problem asks for the speed 'v' at which the gauge pressure at the bottom becomes zero.
# This happens when the initial static pressure (ρ*g*H) is completely converted
# into dynamic pressure (1/2 * ρ * v^2).

# We set up the equation:
# Initial Static Pressure = Dynamic Pressure
# ρ * g * H = (1/2) * ρ * v^2

# We will now solve this for 'v'.

print("To find the speed 'v' where the pressure at the bottom becomes zero, we set the initial static pressure equal to the dynamic pressure.")
print("The equation is: ρ * g * H = (1/2) * ρ * v^2")

print("\nSubstituting the given values into the equation:")
# As requested, here is the final equation with each number
print(f"{rho} * {g} * {H} = (1/2) * {rho} * v^2")

# Show the calculation steps
initial_pressure = rho * g * H
half_rho_term = 0.5 * rho
print(f"{initial_pressure} = {half_rho_term} * v^2")

v_squared = initial_pressure / half_rho_term
print(f"v^2 = {initial_pressure} / {half_rho_term}")
print(f"v^2 = {v_squared}")

# Calculate the final speed v by taking the square root
v = math.sqrt(v_squared)
print(f"v = sqrt({v_squared})")

# Print the final result, formatted to two decimal places
print(f"\nTherefore, the required flow speed is {v:.2f} m/s.")
<<<14.14>>>