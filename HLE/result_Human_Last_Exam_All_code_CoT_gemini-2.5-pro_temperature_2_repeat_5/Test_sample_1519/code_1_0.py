import math

# Define the constants from the problem
g = 10      # Acceleration due to gravity in m/s^2
H = 10      # Depth of the river in meters
rho = 1000  # Density of water in kg/m^3 (though it cancels out)

# The initial hydrostatic pressure at the bottom is P_initial = rho * g * H
# The dynamic pressure due to flow speed v is P_dynamic = 0.5 * rho * v^2
# We need to find the speed v where the pressure at the bottom becomes zero.
# This happens when the initial hydrostatic pressure equals the dynamic pressure.
#
# Equation:
# rho * g * H = 0.5 * rho * v^2
#
# By simplifying, we get:
# g * H = 0.5 * v^2
# v^2 = 2 * g * H
# v = sqrt(2 * g * H)

print("To find the flow speed 'v' where the pressure at the bottom becomes zero, we set the initial hydrostatic pressure equal to the dynamic pressure.")
print("The equation is: ρ * g * H = 0.5 * ρ * v²")
print(f"{rho} * {g} * {H} = 0.5 * {rho} * v²\n")

print("Solving for 'v', the formula becomes: v = sqrt(2 * g * H)")

# Show the equation with the specific numbers
print(f"v = sqrt(2 * {g} * {H})")

# Calculate the final value
v_squared = 2 * g * H
v = math.sqrt(v_squared)

print(f"v = sqrt({v_squared})")
print(f"v = {v}\n")

print(f"The final flow speed is {v:.2f} m/s.")

<<<14.14>>>