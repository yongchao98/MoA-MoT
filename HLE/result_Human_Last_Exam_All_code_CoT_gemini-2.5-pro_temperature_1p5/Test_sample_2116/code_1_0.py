import math

# The derived formula for the expected maximum earthquake magnitude is:
# E[M_max] = pi / (2 * log(2))

# Define the constants in the equation
pi_val = math.pi
log_of_2 = math.log(2)
denominator_val = 2

# Calculate the final result
expected_magnitude = pi_val / (denominator_val * log_of_2)

# Print the final equation and the values of its components
print("The final equation for the expected maximum magnitude is: pi / (2 * log(2))")
print(f"pi = {pi_val}")
print(f"2 = {denominator_val}")
print(f"log(2) = {log_of_2}")
print("-" * 30)
# Print the final result
print(f"Expected maximum earthquake magnitude: {expected_magnitude}")
