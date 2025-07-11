import math

# Define the constants given in the problem
rho = 1000  # Density of water in kg/m^3
g = 10      # Acceleration due to gravity in m/s^2
H = 10      # Depth of the river in meters

# 1. Calculate the initial static gauge pressure at the bottom of the river
initial_pressure = rho * g * H

# 2. Set up the equation based on Bernoulli's principle
# The initial pressure is converted to dynamic pressure (1/2 * rho * v^2)
# We need to find v when the final pressure is 0.
# So, initial_pressure = (1/2) * rho * v^2

# 3. Solve for v
# v^2 = (2 * initial_pressure) / rho
# v = sqrt((2 * initial_pressure) / rho)
v_squared = (2 * initial_pressure) / rho
v = math.sqrt(v_squared)

# 4. Print the explanation and the final result
print("The problem is to find the flow speed 'v' at which the gauge pressure at the bottom of the river becomes zero.")
print("This occurs when the initial static pressure (P) is completely converted to dynamic pressure (1/2 * ρ * v^2).")
print("\nThe equation to solve is: P = (1/2) * ρ * v^2")
print("\nSubstituting the given values:")
print(f"P = ρ * g * H = {rho} * {g} * {H} = {initial_pressure} N/m^2")
print(f"ρ = {rho} kg/m^3")
print("\nThe equation becomes:")
print(f"{initial_pressure} = 0.5 * {rho} * v^2")
print(f"\nSolving for v:")
print(f"v^2 = 2 * {initial_pressure} / {rho} = {v_squared}")
print(f"v = sqrt({v_squared})")
print(f"The flow speed v is: {v} m/s")

# Return the numerical answer in the specified format
# The calculation is sqrt(200)
final_answer = 14.142135623730951
# The problem asks for the answer directly.
# <<<14.14>>> is a good precision. Let's provide more digits for accuracy.
print(f'<<<{v}>>>')