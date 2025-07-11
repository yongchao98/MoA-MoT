import math

# Define the given parameters
rho = 1000  # Density of water in kg/m^3
g = 10      # Acceleration due to gravity in m/s^2
H = 10      # Depth of the river in meters

# The problem states that the initial pressure at the bottom is P_initial = ρ * g * H.
# When the water flows, the pressure at the bottom, P_final, decreases.
# The relationship is: P_final = P_initial - (1/2) * ρ * v^2
# We need to find the speed 'v' where the pressure P_final becomes 0.

# The equation becomes: 0 = (ρ * g * H) - (1/2) * ρ * v^2
# This can be rearranged to: (1/2) * ρ * v^2 = ρ * g * H

print("The equation relating initial pressure and dynamic pressure is:")
print("(1/2) * ρ * v^2 = ρ * g * H")
print("\nSubstituting the numerical values:")
# Note: In the final equation, rho (1000) appears on both sides and can be cancelled out.
# For clarity as requested, we show the full equation first.
print(f"(1/2) * {rho} * v^2 = {rho} * {g} * {H}")

print("\nSimplifying by dividing both sides by ρ (1000) gives:")
print(f"(1/2) * v^2 = {g} * {H}")

# Solve for v^2
# v^2 = 2 * g * H
v_squared = 2 * g * H
print(f"v^2 = 2 * {g} * {H} = {v_squared}")

# Solve for v
v = math.sqrt(v_squared)
print(f"v = sqrt({v_squared})")

print(f"\nTherefore, the flow speed 'v' at which the pressure at the bottom becomes zero is {v:.2f} m/s.")

print(f"\n<<<{v:.2f}>>>")