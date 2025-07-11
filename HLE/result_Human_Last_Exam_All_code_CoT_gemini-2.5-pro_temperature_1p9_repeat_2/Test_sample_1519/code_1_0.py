import math

# Define the given variables
rho = 1000  # Density of water in kg/m^3
g = 10      # Acceleration due to gravity in m/s^2, as given in the problem
H = 10      # Depth of the river in meters

# The problem asks for the flow speed 'v' at which the pressure at the bottom becomes zero.
# According to Bernoulli's principle, the static pressure in a stationary fluid gets converted
# to dynamic pressure (1/2 * rho * v^2) when the fluid starts to flow.

# The initial static gauge pressure at the bottom is P_static = rho * g * H
P_static = rho * g * H

# For the pressure to drop to zero, the entire static pressure must be converted to dynamic pressure.
# So, we set P_static equal to the dynamic pressure term.
# Equation: rho * g * H = 1/2 * rho * v^2

# We can simplify this by canceling out rho (density) from both sides:
# g * H = 1/2 * v^2

# Now, we solve for v:
# v^2 = 2 * g * H
# v = sqrt(2 * g * H)

# Calculate the value of v
v = math.sqrt(2 * g * H)

print("The relationship derived from Bernoulli's principle is:")
print("Initial Static Pressure = Dynamic Pressure")
print("ρ * g * H = 1/2 * ρ * v²")
print("\nSolving for the flow speed v gives:")
print("v = sqrt(2 * g * H)")
print("\nPlugging in the given values:")
print(f"v = sqrt(2 * {g} * {H})")
print(f"The calculated flow speed is {v:.2f} m/s.")

# The problem asks to output the final equation with each number.
# Final Equation: sqrt(2 * 10 * 10)
print("\nThe final equation with each number is:")
print(f"sqrt({2} * {g} * {H})")

print(f"\nFinal Answer: {v}")
<<<14.142135623730951>>>