import math

# Based on the step-by-step derivation, the problem simplifies to calculating
# the value of phi(n) = exp(2*n - 4 + 2/n).
# The problem is stated for n >= 5. Since a single value is expected, we
# compute the result for the smallest valid integer, n=5.

# The final equation for the exponent is: exponent = a*n + b + c/n.
# The instruction asks to output each number in this final equation.
a = 2
b = -4
c = 2

# We choose n=5 for the calculation.
n = 5

# Printing the numbers of the final equation as requested.
print(a)
print(b)
print(c)

# Calculate the exponent.
exponent = a * n + b + c / n

# Calculate the final value of phi(n).
result = math.exp(exponent)

# Print the final numerical result.
print(result)