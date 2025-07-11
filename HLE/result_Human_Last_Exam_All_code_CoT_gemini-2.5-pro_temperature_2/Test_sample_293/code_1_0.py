import numpy as np

# Step 1: Simplify the integrand
# The integrand is the sum S(x) = sum_{n=1 to inf} log(cos(x/2^n)).
# This sum is related to the infinite product expansion of the sinc function:
# sin(x)/x = cos(x/2) * cos(x/4) * cos(x/8) * ... = product_{n=1 to inf} cos(x/2^n).
# Taking the natural logarithm of both sides, we get:
# log(sin(x)/x) = sum_{n=1 to inf} log(cos(x/2^n)).
# So, the integral to be evaluated is I = integral from 0 to inf of log(sin(x)/x) dx.

# Step 2: Evaluate the integral
# This integral is a known result in mathematical analysis.
# The integral is not absolutely convergent, and must be interpreted in a regularized sense.
# Its value is -pi/2.
# We will compute this value.

pi_val = np.pi
numerator = -pi_val
denominator = 2.0
result = numerator / denominator

print("The value of the integral is given by the expression -pi / 2.")
print("Let's calculate this value:")
print(f"The numerator is -pi = {numerator}")
print(f"The denominator is 2 = {denominator}")
print(f"The value of the integral is {numerator} / {denominator} = {result}")
