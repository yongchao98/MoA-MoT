import numpy as np
from scipy.optimize import root_scalar

def find_alpha_0():
    """
    Solves the algebraic equation for alpha_0.
    The problem of finding alpha_0 such that F(alpha_0) = 0 reduces to solving
    the equation 2*alpha^4 - 2*alpha^3 - 1 = 0 for its single positive root.
    """
    
    # Define the function whose root we want to find.
    equation_for_alpha = lambda alpha: 2 * alpha**4 - 2 * alpha**3 - 1

    # We can determine a bracket for the root.
    # f(1) = 2 - 2 - 1 = -1
    # f(1.5) = 2*(1.5)**4 - 2*(1.5)**3 - 1 = 10.125 - 6.75 - 1 = 2.375
    # So, the root lies between 1 and 1.5.
    try:
        sol = root_scalar(equation_for_alpha, bracket=[1, 1.5])
        alpha_0 = sol.root
    except (ImportError, ModuleNotFoundError):
        print("SciPy is not installed. Using a simple numerical solver.")
        # Fallback to a simple bisection method if scipy is not available
        a, b = 1.0, 1.5
        tol = 1e-9
        while (b - a) / 2.0 > tol:
            mid = (a + b) / 2.0
            if equation_for_alpha(mid) == 0:
                alpha_0 = mid
                break
            elif equation_for_alpha(a) * equation_for_alpha(mid) < 0:
                b = mid
            else:
                a = mid
        alpha_0 = (a + b) / 2.0


    print(f"The problem reduces to finding the positive root of the equation:")
    # Outputting the numbers in the equation
    c4, c3, c1, c0 = 2, -2, 0, -1
    print(f"{c4}*\u03B1\u2074 + ({c3})*\u03B1\u00B3 + ({c0}) = 0")
    
    print(f"\nThe largest (and only positive) value \u03B1\u2080 is: {alpha_0}")
    
    print("\nVerification:")
    result = 2 * alpha_0**4 - 2 * alpha_0**3 - 1
    print(f"{c4} * ({alpha_0:.8f})\u2074 + ({c3}) * ({alpha_0:.8f})\u00B3 + ({c0}) = {result:.2e}")

if __name__ == '__main__':
    find_alpha_0()
