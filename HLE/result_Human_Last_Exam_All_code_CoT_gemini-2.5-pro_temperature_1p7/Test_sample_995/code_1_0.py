import numpy as np
from scipy import integrate
import math

# Step 1: Define the optimal F based on the analytical derivation.
# dP/dF = 0 leads to F_opt = 2 / (1 + (2 + pi/4)^2)
pi = np.pi
F_opt = 2 / (1 + (2 + pi / 4)**2)

# Step 2: Define the integrand for the probability calculation.
# This is P(A wins | r, F) / (1/pi) = arccos(sqrt(2*F*r - F^2) / r)
# A's winning probability, given r and F (for F <= 2r), is (1/pi) * integrand.
def integrand(r, F):
    """
    The function inside the integral, which is arccos(...).
    Note: The overall probability includes a (1/pi) factor.
    """
    # The argument of arccos, corresponding to A's optimal strategy for D.
    arg = np.sqrt(2 * F * r - F**2) / r
    # Handle potential floating point inaccuracies where arg might be slightly > 1
    if arg > 1.0:
        arg = 1.0
    return np.arccos(arg)

# Step 3: Calculate the total probability P(F) for F = F_opt.
# The probability is split into two parts based on the relationship between F and r.
# P(F) = Integral from 0 to F/2 of (1) dr + Integral from F/2 to 1 of P(A wins | r, F) dr
# P(F) = F/2 + (1/pi) * Integral from F/2 to 1 of arccos(...) dr

# The first part of the probability: when r < F/2, P(A wins | r) = 1.
# The integral is from 0 to F_opt/2, which evaluates to F_opt/2.
prob_part1 = F_opt / 2

# The second part requires numerical integration.
# We integrate from F_opt/2 to 1. The integrand function takes F_opt as an argument.
integral_result, integral_error = integrate.quad(integrand, F_opt / 2, 1, args=(F_opt,))

# This integral part is divided by pi to get the probability contribution.
prob_part2 = integral_result / pi

# The total minimized probability of A winning.
P_min = prob_part1 + prob_part2

# Step 4: Calculate the final requested value.
final_value = 1 / P_min
floor_value = math.floor(final_value)

# Output the components of the final calculation as requested.
print(f"The optimal value for F is: F_opt = {F_opt}")
print(f"The final probability is the sum of two terms: P_min = P_1 + P_2")
print(f"Term 1 (from r in [0, F_opt/2]): P_1 = {prob_part1}")
print(f"Term 2 (from r in [F_opt/2, 1]): P_2 = {prob_part2}")
print(f"The minimized probability of A winning is: P(A wins) = {P_min}")
print(f"The reciprocal of this probability is: 1 / P(A wins) = {final_value}")
print(f"The floor of the reciprocal is: floor(1 / P(A wins)) = {floor_value}")

print(f"\nFinal equation with numbers:")
print(f"P(A wins) = {prob_part1} + {prob_part2} = {P_min}")
print(f"Result = floor(1 / {P_min}) = floor({final_value}) = {floor_value}")
<<<2>>>