import math

# The supremum of X is given by the expression 8 / (4 * pi^2 + 5)
# This script calculates the numerical value of this expression.

# Value of pi
pi = math.pi

# Calculate pi squared
pi_squared = pi**2

# Calculate 4 * pi^2
four_pi_squared = 4 * pi_squared

# Calculate the denominator 4 * pi^2 + 5
denominator = four_pi_squared + 5

# Calculate the final result
result = 8 / denominator

print("The formula for the supremum of X is: 8 / (4 * pi^2 + 5)")
print(f"To calculate this, we use the following values:")
print(f"pi = {pi}")
print(f"pi^2 = {pi_squared}")
print(f"4 * pi^2 = {four_pi_squared}")
print(f"The denominator (4 * pi^2 + 5) is: {denominator}")
print(f"The supremum of X = 8 / {denominator}")
print(f"Final Answer: {result}")
