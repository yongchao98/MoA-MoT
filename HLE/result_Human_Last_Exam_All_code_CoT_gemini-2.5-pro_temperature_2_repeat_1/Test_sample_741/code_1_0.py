import numpy as np
from scipy.special import iv
from scipy.optimize import root_scalar

def solve_for_x():
    """
    This function finds the largest value of x for which the summation converges to 0.
    The summation is equivalent to the modified Bessel function I_{x-1}(2).
    We need to find the largest root of I_v(2) = 0, where v = x - 1.
    """
    
    # Define the function whose root we want to find.
    # f(v) = I_v(2)
    bessel_func = lambda v: iv(v, 2.0)
    
    # The roots v are all less than -1. We are looking for the largest one.
    # By plotting or numerical evaluation, we find the largest root is between -2.5 and -2.0.
    # iv(-2.5, 2) is negative, while iv(-2.0, 2) is positive.
    bracket = [-2.5, -2.0]
    
    # Use a numerical solver to find the root v.
    solution = root_scalar(bessel_func, bracket=bracket, method='brentq')
    v_root = solution.root
    
    # Calculate x from v using x = v + 1.
    x_root = v_root + 1
    
    # Per instructions, output the numbers in the final equation.
    # The final equation is x = v + 1.
    print(f"The equation to solve is I(v, 2) = 0, where v = x - 1.")
    print(f"The largest root found for v is: {v_root}")
    print(f"The final calculation for x is: {v_root} + 1 = {x_root}")

solve_for_x()