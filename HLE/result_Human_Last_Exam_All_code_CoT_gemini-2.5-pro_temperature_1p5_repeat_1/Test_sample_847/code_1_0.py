import numpy as np
from scipy.optimize import brentq

def equation_for_p_star(p):
    """
    This function represents the condition f'(p)=0, where f(p) is the cost-per-bit.
    Its root is the optimal probability p* that minimizes the cost per bit.
    The equation is derived using calculus on f(p) = (3-2p) / H(p).
    We use natural logarithms (ln) for the calculation. The root is independent of the log base.
    """
    # Using a try-except block to handle p=0 or p=1 which cause log(0) errors.
    try:
        # Prevent log(0) errors at the boundaries of the search interval.
        if p <= 0 or p >= 1:
            return float('inf')
        # The equation from d/dp(f(p))=0 simplifies to the following:
        return 2 * (p * np.log(p) + (1 - p) * np.log(1 - p)) - (3 - 2 * p) * (np.log(1 - p) - np.log(p))
    except ValueError:
        return float('inf')

# We find the optimal probability p* by finding the root of the equation.
# Analysis shows the root lies in the interval (0.5, 1.0). A tighter bound is [0.6, 0.7].
p_star = brentq(equation_for_p_star, 0.6, 0.7)

# The minimal cost per bit, C, can be expressed by a simplified formula
# derived from the f'(p)=0 condition: C = 2*ln(2) / ln(p*/(1-p*)).
# Now, we calculate the numbers that make up this final equation.
numerator = 2 * np.log(2)
denominator = np.log(p_star) - np.log(1 - p_star)
C_min = numerator / denominator

print("The minimal cost to sort n items is asymptotically C * n * log2(n).")
print("The constant factor C is the minimal cost per bit of information.")
print("The final equation for C is: C = (2 * ln(2)) / ln(p* / (1-p*))")
print(f"Optimal probability, p* = {p_star:.5f}")
print(f"Numerator, 2 * ln(2) = {numerator:.5f}")
print(f"Denominator, ln(p* / (1-p*)) = {denominator:.5f}")
print(f"The minimal cost factor, C = {C_min:.3f}")
