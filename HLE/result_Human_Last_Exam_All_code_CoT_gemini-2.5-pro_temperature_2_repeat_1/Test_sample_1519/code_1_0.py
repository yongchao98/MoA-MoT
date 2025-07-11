import math

# 1. Define the given constants for the problem.
# Density of water (rho), though it will cancel out.
rho = 1000  # kg/m^3
# Acceleration due to gravity (g).
g = 10      # m/s^2
# Depth of the river (H).
H = 10      # m

# 2. Set up the equation based on Bernoulli's principle.
# For the pressure at the bottom to become zero, the initial static pressure
# must be completely converted into dynamic pressure.
# Static Pressure = Dynamic Pressure
# rho * g * H   = 1/2 * rho * v^2
# The density 'rho' can be canceled from both sides, which simplifies the equation.
# g * H = 1/2 * v^2

# 3. Solve for the velocity 'v'.
# v^2 = 2 * g * H
# v = sqrt(2 * g * H)

# 4. Perform the calculation and print the steps and result.
print("Based on Bernoulli's principle, for the pressure at the bottom to drop to zero,")
print("the static pressure must equal the dynamic pressure.")
print("\nThe equation is: rho * g * H = 1/2 * rho * v^2")
print("This simplifies to: g * H = 1/2 * v^2")

print("\nNow, we substitute the values into the equation to solve for v:")
# The following print statements will show the equation with the numbers.
print(f"{g} * {H} = (1/2) * v^2")

v_squared = 2 * g * H
print(f"{g * H} = 0.5 * v^2")
print(f"v^2 = 2 * {g * H}")
print(f"v^2 = {v_squared}")

# Calculate the final velocity v
v = math.sqrt(v_squared)
print(f"v = sqrt({v_squared})")
print(f"\nThe required flow speed is {v} m/s.")