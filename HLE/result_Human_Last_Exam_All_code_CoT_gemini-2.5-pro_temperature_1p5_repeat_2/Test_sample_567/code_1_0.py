import math

# Define the integer numbers that form the expression for 'a'
# The expression is a = (7 + 3 * sqrt(5)) / 2
n1 = 7
n2 = 3
n3 = 5
n4 = 2

# Calculate the value of 'a'
sqrt_of_5 = math.sqrt(n3)
numerator = n1 + n2 * sqrt_of_5
a = numerator / n4

# Print the final equation with the numbers, as requested
print(f"The value 'a' where the volume constraint becomes the only obstruction is given by the equation:")
print(f"a = ({n1} + {n2} * sqrt({n3})) / {n4}")

# Print the calculated numerical value
print("\nThe calculated value of 'a' is:")
print(a)