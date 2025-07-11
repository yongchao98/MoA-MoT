import numpy as np
from scipy.integrate import quad

def F_D(d, x1):
    """
    Calculates the Cumulative Distribution Function (CDF) of the distance D = |X - x1|,
    where X is a U[0, 1] random variable.
    d: The distance.
    x1: The fixed point.
    """
    if d < 0:
        return 0
    # The CDF is the length of the interval [x1-d, x1+d] intersected with [0, 1].
    return min(1, x1 + d) - max(0, x1 - d)

def prob_Y_is_X2(x1, x2):
    """
    Calculates the probability that Y = X2 given X1=x1 and X2=x2.
    This happens if |X2-x1| is the median of the three distances from X1.
    """
    d2 = abs(x2 - x1)
    # The probability is 2 * F(d) * (1 - F(d)), where F is the CDF of the other distances.
    f_d = F_D(d2, x1)
    return 2 * f_d * (1 - f_d)

def E_Y2_integrand(x2):
    """
    The outer integrand for calculating E[Y^2].
    This integrates the probability over all possible values of X1.
    """
    # The inner integral is over x1
    inner_integral, _ = quad(lambda x1: prob_Y_is_X2(x1, x2), 0, 1)
    # The full integrand for E[Y^2] calculation
    return 3 * x2**2 * inner_integral

def solve_variance():
    """
    Calculates the variance of Y by numerically integrating to find E[Y^2].
    """
    # E[Y] is 1/2 by symmetry
    e_y = 0.5

    # Numerically integrate to find E[Y^2]
    # E[Y^2] = 3 * E[X2^2 * 1_{Y=X2}] = integral from 0 to 1 of E_Y2_integrand(x2) dx2
    e_y2, err = quad(E_Y2_integrand, 0, 1, limit=100)

    if err > 1e-6:
        print(f"Warning: Numerical integration error ({err}) might be high.")

    # Variance(Y) = E[Y^2] - (E[Y])^2
    var_y = e_y2 - e_y**2

    print("Determining the variance of Y, Var(Y) = E[Y^2] - (E[Y])^2")
    print("-" * 50)
    print(f"By symmetry, the expected value E[Y] = {e_y}")
    print(f"The expected value of Y^2, E[Y^2], is computed numerically: {e_y2:.8f}")
    print(f"The variance is Var(Y) = {e_y2:.8f} - ({e_y})^2")
    print(f"                     = {e_y2:.8f} - {e_y**2:.8f}")
    print(f"                     = {var_y:.8f}")
    print("\nThe exact value of the variance is 1/24.")
    print(f"1/24 is approximately {1/24:.8f}")

solve_variance()