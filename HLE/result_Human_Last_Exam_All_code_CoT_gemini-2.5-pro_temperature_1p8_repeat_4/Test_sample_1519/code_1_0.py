import math

# Description of the problem and variables
print("This script calculates the flow speed 'v' at which the pressure at the bottom of a river becomes zero.")
print("The principle used is the conversion of hydrostatic pressure to dynamic pressure.\n")

# Given values
rho = 1000  # density of water in kg/m^3
g = 10      # acceleration due to gravity in m/s^2
H = 10      # depth of the river in meters

# The hydrostatic gauge pressure at the bottom of the static river
static_pressure = rho * g * H

print("First, we state the equation equating hydrostatic and dynamic pressure:")
print("rho * g * H = (1/2) * rho * v^2")
print("Plugging in the known values for rho, g, and H:")
print(f"{rho} * {g} * {H} = (1/2) * {rho} * v^2")
print("The left side of the equation calculates the initial static pressure:")
print(f"{static_pressure} = 0.5 * {rho} * v^2\n")

print("Now, we solve for v^2:")
print(f"v^2 = (2 * {static_pressure}) / {rho}")
v_squared = (2 * static_pressure) / rho
print(f"v^2 = {v_squared}\n")

print("Finally, we solve for v by taking the square root:")
v = math.sqrt(v_squared)
print(f"v = sqrt({v_squared})")
print(f"v = {v}")

print("\nTherefore, the flow speed at which the pressure at the bottom decreases to zero is approximately 14.14 m/s.")