# This script defines the coefficients and exponents for the formula of P(n)
# and uses them to construct and print the formula string, outputting each number.

# The first term's components in the exponent come from the B_4 Bernoulli term.
num1 = 1
den1 = 720
power1 = 2

# The second term's components in the exponent come from the B_6 Bernoulli term.
num2 = 1
den2 = 5040
power2 = 4

# Construct and print the formula for P(n).
print(f"P(n) = exp({num1}/({den1}*n^{power1}) - {num2}/({den2}*n^{power2}))")