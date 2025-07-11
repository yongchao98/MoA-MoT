import numpy as np
from scipy.optimize import brentq

# The problem of finding the minimal cost per bit of information boils down to
# solving a cubic equation for the optimal probability 'p' of a 'yes' answer.

# 1. Define and solve the equation for the optimal probability 'p'.
# The equation, derived from balancing worst-case costs, is p^3 + p - 1 = 0.
def p_equation(p):
    # The coefficients are 1, 1, -1.
    return p**3 + p - 1

# We find the real root of this equation in the interval (0, 1).
p_optimal = brentq(p_equation, 0, 1)

# 2. Calculate the minimal cost 'C' per bit of information.
# The cost constant C is given by the formula C = -1 / log2(p).
cost_per_bit = -1 / np.log2(p_optimal)

# 3. Output the results, showing the numbers from the equations as requested.
print("The optimal probability 'p' is the real root of the equation: (1)*p^3 + (1)*p + (-1) = 0.")
print(f"The numerical value for p is: {p_optimal:.5f}")

print("\nThe minimal cost 'C' per bit of information is calculated using the formula: C = (-1) / log2(p).")
print(f"The resulting minimal cost factor, rounded to 3 decimal places, is: {cost_per_bit:.3f}")
