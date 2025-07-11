import numpy as np

# Based on the analysis, the exact value of l(a) is (a-1) * (ln(2) + 1).
# The problem asks to output each number in the final equation.
# The equation is l(a) = (a - 1) * (ln(2) + 1).
# The numbers are 1, 2, and 1.

a_minus_val = 1
log_argument = 2
plus_val = 1

# We construct and print the final equation string.
final_equation = f"l(a) = (a - {a_minus_val}) * (ln({log_argument}) + {plus_val})"
print(final_equation)

# We can also print the numerical value of the constant factor.
constant_factor = np.log(log_argument) + plus_val
print(f"l(a) = (a - {a_minus_val}) * {constant_factor}")
