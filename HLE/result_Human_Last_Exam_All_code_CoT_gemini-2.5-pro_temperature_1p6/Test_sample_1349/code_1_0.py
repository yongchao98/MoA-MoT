import numpy as np

# The supremum of X is given by the expression 8 / (4*pi^2 + 5).
# Here are the numbers in the final equation:
num_8 = 8
num_4 = 4
num_5 = 5
pi = np.pi

# Calculate the value of the expression
sup_X = num_8 / (num_4 * pi**2 + num_5)

# Output the equation and the result
print(f"The final expression for the supremum of X is {num_8} / ({num_4}*pi^2 + {num_5})")
print(f"The numerical value is: {sup_X}")
