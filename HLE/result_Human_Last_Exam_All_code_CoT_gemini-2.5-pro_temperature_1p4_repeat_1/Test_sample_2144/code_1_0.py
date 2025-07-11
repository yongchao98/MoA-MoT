import numpy as np

def solve_trajectory():
    """
    This script solves for the position x0 where y(x0)=-3.
    The problem is solved by parameterizing the trajectory by its slope p=dy/dx.
    This leads to a quartic equation for p when y=-3.
    """
    
    # The condition y = -3 leads to the quartic equation for p: p^4 - 18p - 27 = 0.
    # We define the coefficients of the polynomial p^4 + 0p^3 + 0p^2 - 18p - 27.
    coeffs = [1, 0, 0, -18, -27]
    
    # Find the roots of the polynomial equation.
    roots = np.roots(coeffs)
    
    print("The trajectory parameter 'p' must satisfy the quartic equation: p^4 - 18p - 27 = 0")
    
    # The parametric solution requires the term sqrt(2p+3), so 2p+3 must be non-negative.
    # We filter for real roots that satisfy this condition (p >= -1.5).
    valid_roots = [r.real for r in roots if np.isreal(r) and 2 * r.real + 3 >= 0]
    
    # The particle starts at y=-1 and moves to y=-3. Analysis of the function y(p)
    # shows that y decreases only when p > 0.
    # Therefore, we must choose the positive real root.
    p_solution = 0
    for r in valid_roots:
        if r > 0:
            p_solution = r
            break

    # The root finding might give a float, so we round to the nearest integer.
    p_solution = round(p_solution)

    print(f"The required value of the slope is p = {p_solution}.")
    
    # Now, we calculate the position x0 using the parametric equation for x(p).
    # x(p) = -3(p+1) / sqrt(2p+3)
    p = p_solution
    numerator = -3 * (p + 1)
    denominator_val = 2 * p + 3
    denominator = np.sqrt(denominator_val)
    x0 = numerator / denominator
    
    print("\nThe position x0 is calculated using the equation: x0 = -3(p + 1) / sqrt(2p + 3)")
    print(f"Substituting p = {p}:")
    print(f"x0 = -3({p} + 1) / sqrt(2*{p} + 3)")
    print(f"x0 = {numerator} / sqrt({denominator_val})")
    print(f"x0 = {numerator} / {denominator}")
    print(f"x0 = {x0}")
    
solve_trajectory()
<<<-4.0>>>