import math

# The problem is to calculate the integral of u(x,y,-y,1) with respect to x from 0 to 1.
# Based on the properties of the given equation and the initial condition,
# we deduce that u(x,y,-y,t) is time-independent.
# So, we can calculate the integral of u(x,y,-y,0) instead.
# The integral is I = integral from 0 to 1 of -3 * (2*e^(2x) + e^x) / (e^(2x) + e^x + 1) dx.

# This integral can be solved analytically.
# Let f(x) = e^(2x) + e^x + 1. The derivative f'(x) = 2*e^(2x) + e^x.
# The integral is -3 * integral(f'(x)/f(x) dx), which is -3 * [ln(f(x))].
# Evaluating from 0 to 1:
# I = -3 * ln(f(1)) - (-3 * ln(f(0)))
# I = 3 * (ln(f(0)) - ln(f(1)))
# I = 3 * (ln(e^0 + e^0 + 1) - ln(e^2 + e + 1))
# I = 3 * (ln(3) - ln(e^2 + e + 1))
# I = 3 * ln(3 / (e^2 + e + 1))

# Now, we calculate this value using Python.

e = math.e

# Numerator of the argument of the logarithm
num = 3

# Denominator of the argument of the logarithm
den = e**2 + e + 1

# Calculate the final result
result = 3 * math.log(num / den)

# As requested, output the numbers in the final equation and the result
print("The analytical solution to the integral is of the form: 3 * ln(3 / (e^2 + e + 1))")
print("\nCalculating the components:")
print(f"e = {e}")
print(f"e^2 = {e**2}")
print(f"Denominator = e^2 + e + 1 = {den}")
print(f"Final equation: 3 * ln({num} / {den})")
print(f"\nResult: {result}")
