import math

# Define the constants given in the problem
rho = 1000  # Density of water in kg/m^3
g = 10      # Acceleration due to gravity in m/s^2
H = 10      # Depth of the river in meters

# The principle is that the initial hydrostatic pressure (P_hydrostatic)
# is converted into dynamic pressure (P_dynamic) due to the flow.
# We need to find the velocity 'v' for which the gauge pressure at the bottom becomes zero.
# This happens when P_hydrostatic = P_dynamic.
# The equation is: rho * g * H = 0.5 * rho * v**2

# We can solve for v^2: v**2 = 2 * g * H
v_squared = 2 * g * H

# Now, calculate the velocity v by taking the square root
v = math.sqrt(v_squared)

# Print the equation with the numerical values substituted in.
# We show the full equation first for clarity, as requested.
print("The problem is to find the velocity 'v' where the initial hydrostatic pressure equals the dynamic pressure.")
print("The equation is: rho * g * H = (1/2) * rho * v^2")
print(f"Substituting the values: {rho} * {g} * {H} = 0.5 * {rho} * v^2")

# Calculate the pressure on the left side
pressure = rho * g * H
print(f"This simplifies to: {pressure} = {0.5 * rho} * v^2")
print(f"Solving for v^2: v^2 = {pressure} / {0.5 * rho}")
print(f"v^2 = {v_squared}")
print(f"Solving for v: v = sqrt({v_squared})")
print(f"The required flow speed v is: {v} m/s")
