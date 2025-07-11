import math

# Description: This script calculates the flow speed 'v' of a river
# at which the pressure at the bottom becomes zero, based on Bernoulli's principle.

# Given parameters
rho = 1000  # density of water in kg/m^3
g = 10      # acceleration due to gravity in m/s^2
H = 10      # depth of the river in meters

# The initial pressure at the bottom is the hydrostatic pressure P = rho * g * H.
# This pressure is converted into kinetic energy K = (1/2) * rho * v^2.
# We set P = K to find the speed 'v' where the pressure becomes 0.
# The equation is: rho * g * H = (1/2) * rho * v^2

print("The relationship derived from Bernoulli's principle for this problem is:")
print("Initial Pressure = Final Kinetic Energy")
print("rho * g * H = (1/2) * rho * v^2")

# Substitute the given values into the equation
initial_pressure = rho * g * H
half_rho = 0.5 * rho

print("\nSubstituting the given values:")
print(f"{rho} * {g} * {H} = (1/2) * {rho} * v^2")
print(f"{initial_pressure} = {half_rho} * v^2")

# Solve for v^2
# Note: The 'rho' term can be cancelled from both sides, simplifying to g * H = 0.5 * v^2
v_squared = (rho * g * H) / (0.5 * rho)
print("\nSolving for v^2:")
print(f"v^2 = {initial_pressure} / {half_rho}")
print(f"v^2 = {v_squared}")

# Solve for v by taking the square root
v = math.sqrt(v_squared)
print("\nSolving for v:")
print(f"v = sqrt({v_squared})")
print(f"v = {v}")

print(f"\nTherefore, the flow speed at which the pressure at the bottom decreases to zero is {v} m/s.")
<<<14.142135623730951>>>