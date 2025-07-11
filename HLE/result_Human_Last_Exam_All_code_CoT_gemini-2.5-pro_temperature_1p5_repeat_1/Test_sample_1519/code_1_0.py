import math

# Define the given physical constants and parameters
rho = 1000  # Density of water in kg/m^3
g = 10      # Acceleration due to gravity in m/s^2
H = 10      # Depth of the river in meters

# Calculate the initial hydrostatic pressure at the bottom of the river
initial_pressure = rho * g * H

print(f"The initial gauge pressure at the bottom of the stationary river is P = rho * g * H:")
print(f"P = {rho} * {g} * {H} = {initial_pressure} N/m^2")
print("-" * 40)

# Explain the setup for finding the flow speed v
print("When the water flows, the pressure at the bottom (P_new) is reduced by the dynamic pressure:")
print("P_new = (rho * g * H) - (1/2 * rho * v^2)")
print("\nWe want to find the speed v where the pressure at the bottom becomes 0.")
print("Setting P_new = 0, we get the equation:")
print(f"0 = {initial_pressure} - (1/2 * {rho} * v^2)")

# Rearranging the equation to solve for v gives: v = sqrt(2 * g * H)
print("\nSolving for v:")
print("v = sqrt(2 * g * H)")
print("-" * 40)

# Calculate the final value for v and print the result with the numbers in the equation
v_squared = 2 * g * H
v = math.sqrt(v_squared)

print("Plugging in the values to find the speed v:")
print(f"v = sqrt(2 * {g} * {H})")
print(f"v = sqrt({v_squared})")
print(f"v = {v} m/s")

# Final answer block
# The question requires the numerical answer in a specific format
final_answer = v
# print(f"\n<<<{final_answer}>>>") # This would be the final output line
# As per instructions, let's output it for the user.
print("\nThe flow speed v at which the pressure at the bottom decreases to zero is:")
print(final_answer)
print("<<<14.142135623730951>>>")