import math

# Define the constants given in the problem
rho = 1000  # Density of water in kg/m^3
g = 10      # Acceleration due to gravity in m/s^2
H = 10      # Depth of the river in meters

# The initial pressure at the bottom of the static river is P_initial = ρ * g * H.
# This pressure is converted to kinetic energy (1/2 * ρ * v^2) when the water flows.
# We want to find the velocity 'v' for which the pressure at the bottom becomes zero.
# The equation is: ρ * g * H = 1/2 * ρ * v^2

# Calculate the initial pressure for clarity in the output
initial_pressure = rho * g * H

print("To find the flow speed 'v' at which the pressure at the bottom becomes zero, we set the initial static pressure equal to the kinetic energy term.")
print("The equation is: ρ * g * H = (1/2) * ρ * v^2")
print("\nPlugging in the values, we get the final equation:")
# We need to output each number in the final equation
print(f"{rho} * {g} * {H} = 0.5 * {rho} * v^2")

# Simplifying the equation:
# g * H = 0.5 * v^2
# v^2 = 2 * g * H
v_squared = 2 * g * H
v = math.sqrt(v_squared)

print(f"\nSolving for v:")
print(f"v^2 = 2 * {g} * {H}")
print(f"v^2 = {v_squared}")
print(f"v = sqrt({v_squared})")
print(f"The flow speed v is {v:.2f} m/s.")

# Final answer format
# The final result is sqrt(200) which is approximately 14.14
final_answer = v
# The prompt requests the answer in a specific format at the end.
# I will output the final numerical answer as requested.
# The value is math.sqrt(200) = 14.1421356...
# Let's round to a reasonable number of decimal places or provide the exact form if simpler.
# 14.14 seems appropriate.