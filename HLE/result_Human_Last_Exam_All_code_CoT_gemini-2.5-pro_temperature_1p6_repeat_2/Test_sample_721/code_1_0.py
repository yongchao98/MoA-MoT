import numpy as np
from scipy.integrate import quad
import math

# For function 2, the inequality holds if (2 - sqrt(2)) * I_1 <= I_2
# where I_1 is the integral for S_1 and I_2 is the integral for |a_0|

# Integrand for S_1 = integral from 0 to 1 of 1/sqrt(1-x^4) dx
def integrand_I1(x):
    return 1.0 / np.sqrt(1 - x**4)

# Integrand for |a_0| = integral from 0 to 1 of 1/sqrt(t+t^3) dt
def integrand_I2(t):
    # The integral is improper at t=0, but integrable.
    # It can be handled by quad, but for robustness let's check.
    if t == 0:
        return np.inf # scipy.quad can often handle this
    return 1.0 / np.sqrt(t + t**3)

# Perform the numerical integration
# For I_1, the upper limit is 1, where the integrand is singular.
I1, err1 = quad(integrand_I1, 0, 1)
# For I_2, the lower limit is 0, where the integrand is singular.
I2, err2 = quad(integrand_I2, 0, 1)

# The condition to check
lhs_value = (2 - math.sqrt(2)) * I1
rhs_value = I2

# Output the plan and the results of the check
print("For function 2, we are checking the condition (2 - sqrt(2)) * I_1 <= I_2")
print(f"where I_1 = integral from 0 to 1 of 1/sqrt(1-x^4) dx")
print(f"and I_2 = |a_0| = integral from 0 to 1 of 1/sqrt(t+t^3) dt")
print("-" * 30)
print(f"Numerically computed I_1: {I1:.4f}")
print(f"Numerically computed I_2 (|a_0|): {I2:.4f}")
print("-" * 30)

print("Now checking the inequality:")
print(f"LHS of condition: (2 - sqrt(2)) * I_1 = (2 - {math.sqrt(2):.4f}) * {I1:.4f} = {lhs_value:.4f}")
print(f"RHS of condition: I_2 = {rhs_value:.4f}")
print("-" * 30)

if lhs_value <= rhs_value:
    print("The condition is satisfied: TRUE.")
    print(f"Therefore, function 2 satisfies the original inequality.")
else:
    print("The condition is not satisfied: FALSE.")
    print(f"Therefore, function 2 does not satisfy the original inequality.")
