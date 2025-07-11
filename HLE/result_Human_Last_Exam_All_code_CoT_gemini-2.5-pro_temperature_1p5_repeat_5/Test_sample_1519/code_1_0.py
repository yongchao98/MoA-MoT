import math

# Define the constants given in the problem
g = 10  # Acceleration due to gravity in m/s^2
H = 10  # Depth of the river in meters

# According to the problem statement, we want to find the speed 'v'
# where the pressure at the bottom of the river becomes 0.
# This occurs when the initial static pressure (P_static) equals the
# dynamic pressure (P_dynamic) caused by the flow.

# P_static = ρ * g * H
# P_dynamic = (1/2) * ρ * v^2

# Setting them equal:
# ρ * g * H = (1/2) * ρ * v^2

# The density 'ρ' cancels out, simplifying the equation to:
# g * H = (1/2) * v^2

# We can now solve for v:
# v^2 = 2 * g * H
# v = sqrt(2 * g * H)

# Now, we substitute the values of g and H into the equation to find v.
v_squared = 2 * g * H
v = math.sqrt(v_squared)

print("To find the flow speed v, we set the initial static pressure equal to the dynamic pressure:")
print("ρ * g * H = (1/2) * ρ * v^2")
print("Solving for v gives the equation: v = sqrt(2 * g * H)")
print("\nSubstituting the given values into the equation:")
print(f"v = sqrt(2 * {g} * {H})")
print(f"v = sqrt({v_squared})")
print(f"The flow speed v is {v:.2f} m/s.")
