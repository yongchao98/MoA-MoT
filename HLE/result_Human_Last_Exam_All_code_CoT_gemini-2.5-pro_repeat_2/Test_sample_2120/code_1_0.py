# The locations of the poles are given by the roots of the denominator
# of the expression for the B-field.
# The polynomial is 4*z^4 - z^3 + z^2 + 1 = 0.

# We use Vieta's formulas to find the sum of the roots without solving the equation.
# For a polynomial a_n*z^n + a_{n-1}*z^(n-1) + ... + a_0 = 0, the sum of roots
# is -a_{n-1}/a_n.

# Identify coefficients from the polynomial: 4z^4 - 1z^3 + 1z^2 + 1 = 0
a_n = 4  # Coefficient of the highest power term, z^4
a_n_minus_1 = -1  # Coefficient of the second-highest power term, z^3
n = 4  # The degree of the polynomial, which is the number of roots.

# Calculate the sum of the roots
sum_of_roots = -a_n_minus_1 / a_n

# Calculate the average value of the roots
average_value = sum_of_roots / n

# The final equation for the average is (sum of roots) / n.
# We print the calculation with the specific numbers.
print(f"The equation for the average value is: (-({a_n_minus_1}) / {a_n}) / {n}")
print(f"The calculated average value is: {average_value}")