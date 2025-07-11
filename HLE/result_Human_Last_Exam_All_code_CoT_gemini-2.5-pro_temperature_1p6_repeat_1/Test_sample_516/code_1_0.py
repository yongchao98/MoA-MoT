import math

# The expected value of T is given by the equation E[T] = a - b * e^(c/d)
# Here are the numbers from the final equation:
a = 3
b = 2
c = 1
d = 3

# Calculate the final result
result = a - b * math.exp(c / d)

# Print the final equation and its value
print(f"The derived equation for the expected value is: E[T] = {a} - {b} * e^({c}/{d})")
print(f"The calculated value is: {result}")