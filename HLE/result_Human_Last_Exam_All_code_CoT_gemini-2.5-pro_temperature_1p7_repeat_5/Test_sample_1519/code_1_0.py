import math

# Define the given physical constants
rho = 1000  # Density of water in kg/m^3
g = 10      # Acceleration due to gravity in m/s^2
H = 10      # Depth of the river in m

# The problem asks for the speed 'v' at which the pressure at the bottom becomes zero.
# This occurs when the initial hydrostatic pressure (rho * g * H) is completely
# converted into dynamic pressure (0.5 * rho * v^2).

# We set up the equation:
# rho * g * H = 0.5 * rho * v^2

# We solve this equation for v:
# v = sqrt(2 * g * H)
v = math.sqrt(2 * g * H)

# Print the final equation with the numerical values substituted
print("The relationship between the initial pressure and the dynamic pressure is:")
print("rho * g * H = 0.5 * rho * v^2")
print("\nSubstituting the values into the equation:")
print(f"{rho} * {g} * {H} = 0.5 * {rho} * v^2")

# Print the result
print(f"\nSolving for v, the required flow speed is {v:.2f} m/s.")

<<<14.14>>>