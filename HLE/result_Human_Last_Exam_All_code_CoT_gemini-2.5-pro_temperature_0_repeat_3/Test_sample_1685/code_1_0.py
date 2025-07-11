import math

# The asymptotic solution is of the form y(x) = A / x.
# From the dominant balance analysis, the coefficient A is -cbrt(6).

# Calculate the value of the coefficient A
A = -(6**(1/3))

# Round the coefficient to two decimal places as requested
A_rounded = round(A, 2)

# The exponent of x is -1
n = -1

# The problem asks to output each number in the final equation.
# We will print the final expression in a clear format.
# The expression is y(x) = A_rounded / x, which is A_rounded * x**n
print(f"The approximate solution for large x is y(x) = {A_rounded} / x")
print("The numbers in this equation are:")
print(f"Coefficient A: {A_rounded}")
print(f"Exponent n: {n}")
