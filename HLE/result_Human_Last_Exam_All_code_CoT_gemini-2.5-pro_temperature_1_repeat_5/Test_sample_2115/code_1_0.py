import math

# The problem simplifies to the evaluation of a definite integral, which is solved analytically.
# The final expression for the integral is 3 * log(3 / (e^2 + e + 1)).
# This script calculates the numerical value of this expression.

# The numbers in the final equation are 3, e, 2, and 1.
# 'e' is the mathematical constant, math.e in Python.
final_result = 3 * math.log(3 / (math.e**2 + math.e + 1))

# Print the final equation with numerical values and the result.
print(f"3 * log(3 / ({math.e}**2 + {math.e} + 1)) = {final_result}")