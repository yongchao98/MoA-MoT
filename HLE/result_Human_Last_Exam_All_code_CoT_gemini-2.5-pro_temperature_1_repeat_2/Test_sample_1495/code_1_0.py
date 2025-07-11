import math

# This problem is a classic in geometric probability. The probability (P) is the
# expected value of the ratio of the inner triangle's area to the outer
# triangle's area. The analytical solution for this is a known result.
# The formula for the probability is: P = 10 - π^2

# Define the numbers used in the final equation.
number_ten = 10
pi_value = math.pi

# Calculate pi squared.
pi_squared = pi_value ** 2

# Calculate the final probability.
probability = number_ten - pi_squared

# The user prompt asks to output each number in the final equation.
# The equation is: Probability = 10 - π^2
print("The final equation for the probability (P) is:")
print(f"P = {number_ten} - (π)^2")

print("\nBreaking down the calculation:")
print(f"The first number is: {number_ten}")
print(f"The value of π is: {pi_value}")
print(f"The value of π squared is: {pi_squared}")

print("\nSubstituting the values into the equation:")
print(f"P = {number_ten} - {pi_squared}")
print(f"P = {probability}")
