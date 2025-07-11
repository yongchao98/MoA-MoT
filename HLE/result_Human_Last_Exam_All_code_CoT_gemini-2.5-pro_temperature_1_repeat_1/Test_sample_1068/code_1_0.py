import math

# Plan:
# The step-by-step derivation shows that the complex function l(p) simplifies to the expression 2 * p * log(2).
# The problem asks for the value of l(14).
# This script will substitute p=14 into the simplified formula and print the result.
# It will also display the components of the final calculation as requested by the prompt.

# The value of p as given in the problem
p = 14

# The simplified formula for l(p) is 2 * p * log(2).
# Calculate the coefficient for log(2)
coefficient = 2 * p

# Calculate the natural logarithm of 2
log_of_2 = math.log(2)

# Compute the final result
result = coefficient * log_of_2

# Output the explanation and the final equation with its components
print("The complex problem simplifies to the formula l(p) = 2 * p * log(2).")
print(f"For the given value p = {p}, the calculation is:")
print(f"{coefficient} * log(2) = {coefficient} * {log_of_2}")
print(f"Final value = {result}")
