import math

# Define the given physical constants
rho = 1000  # Density of water in kg/m^3
g = 10      # Acceleration due to gravity in m/s^2
H = 10      # Depth of the river in meters

print("This script calculates the flow speed 'v' at which the pressure at the bottom of a river becomes zero.\n")

# Step 1: State the initial pressure at the bottom (hydrostatic pressure).
P_static = rho * g * H
print(f"First, we calculate the initial hydrostatic pressure at the bottom when the water is at rest:")
print(f"P_static = ρ * g * H")
print(f"P_static = {rho} * {g} * {H} = {P_static} N/m^2\n")

# Step 2: Apply Bernoulli's principle.
# For the pressure at the bottom to decrease to zero, the initial static pressure
# must be completely converted into dynamic pressure, (1/2) * ρ * v^2.
print("According to Bernoulli's principle, this static pressure is converted into dynamic pressure as the water flows.")
print("We set the static pressure equal to the dynamic pressure to find the required speed 'v':")
print("P_static = (1/2) * ρ * v^2")
print(f"ρ * g * H = (1/2) * ρ * v^2\n")


# Step 3: Solve the equation for v.
# We can simplify the equation by canceling out the density ρ.
# g * H = (1/2) * v^2
# v^2 = 2 * g * H
# v = sqrt(2 * g * H)
print("Solving the equation for v:")
print("v = sqrt(2 * g * H)")

v_squared = 2 * g * H
v = math.sqrt(v_squared)

print(f"v = sqrt(2 * {g} * {H})")
print(f"v = sqrt({v_squared})")
print(f"v = {v} m/s\n")

print(f"The flow speed required for the pressure at the bottom to decrease to zero is {v} m/s.")

print(f'<<<{v}>>>')