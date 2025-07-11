import math

# The problem is to evaluate the sum S = sum_{n=2 to inf} (-1)^n * S_n / n,
# where S_n is the n-th harmonic number.
# Through mathematical analysis, the sum can be shown to be exactly equal to the expression:
# S = 1 - (pi^2 / 12) + ((ln(2))^2 / 2)

# Here are the constants from the final equation
c1 = 1
c2 = 12
c3 = 2

# We use the math library for the values of pi and the natural logarithm.
pi = math.pi
ln2 = math.log(2)

# Calculate the individual terms in the expression
term1 = float(c1)
term2 = pi**2 / c2
term3 = ln2**2 / c3

# The final equation is S = term1 - term2 + term3
final_value = term1 - term2 + term3

# As per the instruction, we print the components of the final equation and the final result.
print(f"The final expression is of the form: {c1} - (pi^2 / {c2}) + (ln(2)^2 / {c3})")
print(f"Value of pi^2 / {c2}: {term2}")
print(f"Value of (ln(2))^2 / {c3}: {term3}")
print(f"The final value of the sum is: {final_value}")