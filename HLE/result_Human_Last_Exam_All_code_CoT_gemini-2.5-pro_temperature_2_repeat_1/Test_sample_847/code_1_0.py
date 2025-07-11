import numpy as np
from scipy.optimize import root_scalar

def equation(p):
    """
    Defines the cubic equation p^3 + p - 1 = 0, which determines
    the optimal probability 'p' for a 'yes' answer.
    """
    return p**3 + p - 1

# We solve for the real root of the equation in the interval [0, 1],
# since 'p' represents a probability.
solution = root_scalar(equation, bracket=[0, 1])
p_optimal = solution.root

# The minimal cost per bit of information is given by the formula C = -1 / log2(p).
cost_per_bit = -1 / np.log2(p_optimal)

# The final equation gives the value for the minimal cost factor 'C'.
# We display the numbers that constitute this calculation, as requested.
print(f"Minimal cost factor C = -1 / log2({p_optimal:.4f}) = {cost_per_bit:.3f}")