import math

# Define the constants based on the problem statement
rho = 1000  # Density of water in kg/m^3
g = 10      # Acceleration due to gravity in m/s^2
H = 10      # Depth of the river in meters

# The pressure at the bottom when the water is at rest is P_initial = ρ * g * H
# When the water flows, this pressure is reduced by the dynamic pressure term 1/2 * ρ * v^2.
# The new pressure is P_final = (ρ * g * H) - (1/2 * ρ * v^2).

# We want to find the speed 'v' for which the pressure at the bottom becomes zero.
# Set P_final = 0:
# 0 = (ρ * g * H) - (1/2 * ρ * v^2)
# ρ * g * H = 1/2 * ρ * v^2

# We can simplify by dividing by ρ:
# g * H = 1/2 * v^2

# Now, we solve for v:
# v^2 = 2 * g * H
# v = sqrt(2 * g * H)

# Perform the calculation
v_squared = 2 * g * H
v = math.sqrt(v_squared)

# Print the step-by-step solution as requested
print("To find the flow speed 'v' where the pressure at the bottom becomes zero, we set the initial hydrostatic pressure equal to the dynamic pressure.")
print("The equation is: ρ * g * H = 1/2 * ρ * v^2")
print("This simplifies to v = sqrt(2 * g * H).")
print("\nSubstituting the given values into the equation:")
print(f"v = sqrt(2 * {g} * {H})")
print(f"v = sqrt({v_squared})")
print(f"\nThe required flow speed is {v} m/s.")

print("<<<14.142135623730951>>>")