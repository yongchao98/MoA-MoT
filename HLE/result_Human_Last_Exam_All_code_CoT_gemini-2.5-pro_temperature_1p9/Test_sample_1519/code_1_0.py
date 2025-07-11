import math

# Define the constants given in the problem
g = 10  # acceleration due to gravity in m/s^2
H = 10  # depth of the river in meters

# The equation to solve is v = sqrt(2 * g * H)
# We are finding the speed v where the dynamic pressure (1/2 * rho * v^2)
# equals the initial hydrostatic pressure (rho * g * H).

# Calculate the value of v
v_squared = 2 * g * H
v = math.sqrt(v_squared)

# Print the explanation and the final equation with the numbers plugged in
print("To find the flow speed 'v' at which the pressure at the bottom becomes zero, we set the initial hydrostatic pressure equal to the dynamic pressure.")
print("The equation is: ρ * g * H = 1/2 * ρ * v²")
print("This simplifies to: v = sqrt(2 * g * H)")
print("\nPlugging in the values:")
print(f"v = sqrt(2 * {g} * {H})")
print(f"v = sqrt({v_squared})")
print(f"v = {v:.2f} m/s")

# Final answer format as requested by the user prompt
# The question in the prompt is: At what flow speed v does the pressure at the bottom decrease to zero 0?
# The value is v, which we just calculated. I'll provide the numeric value rounded to 2 decimal places.
# final_answer = f"<<<{v:.2f}>>>"
# The instruction is "directly return the answer", let me check if that implies printing the value inside the angle brackets. It seems so. Let me just provide the numeric value rounded.