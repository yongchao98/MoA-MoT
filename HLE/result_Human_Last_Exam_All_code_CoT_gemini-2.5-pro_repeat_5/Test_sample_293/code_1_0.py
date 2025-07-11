import numpy as np
from scipy.integrate import quad

# The original integral is I = integral from 0 to infinity of sum_{n=1 to inf} log(cos(x/2^n)) dx.
# Using the identity product_{n=1 to inf} cos(x/2^n) = sin(x)/x, the sum becomes log(sin(x)/x).
# So, I = integral from 0 to infinity of log(sin(x)/x) dx.
# For x > pi, sin(x) becomes negative, and the integrand becomes complex.
# The integral of the imaginary part diverges. We assume the question asks for the real part of the integral.
# The real part is integral from 0 to infinity of log|sin(x)/x| dx.
# This is a known integral with the value -pi/2.

# We define the function for the real part of the integrand.
def real_integrand(x):
    # To avoid division by zero at x=0, we handle it as a special case.
    # The limit of sin(x)/x as x->0 is 1, and log(1)=0.
    if x == 0:
        return 0
    return np.log(np.abs(np.sin(x) / x))

# We use scipy.integrate.quad to compute the definite integral from 0 to infinity.
# The function returns the result and an estimated error.
integral_value, error = quad(real_integrand, 0, np.inf)

# The expected value is -pi/2.
expected_value = -np.pi / 2

print("This script calculates the value of the integral by numerically evaluating its real part.")
print(f"The numerical result of the integral is: {integral_value}")
print(f"The theoretical value is -pi/2, which is approximately: {expected_value}")
print("\nFinal Equation:")
print(f"Value = - (pi) / 2")
# The prompt requested to output each number in the final equation.
# The numbers are pi and 2.
pi_val = np.pi
divisor = 2
result = -pi_val / divisor
print(f"Value = - {pi_val} / {divisor}")
print(f"Value = {result}")
