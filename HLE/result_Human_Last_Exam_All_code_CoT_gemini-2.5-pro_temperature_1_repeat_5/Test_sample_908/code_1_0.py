import sympy

def solve_charges():
    """
    Calculates and prints the expressions for total volume and surface charges
    based on the formulas provided in Option B of the problem.
    """
    # Define symbolic variables
    V, epsilon, pi, L, a, b = sympy.symbols('V varepsilon pi L a b', real=True, positive=True)

    # Expressions from Option B
    q_v_expr = -4 * V * epsilon * pi * L / (1 - a**2 / b**2)
    q_s_a_expr = 2 * pi * L * V * epsilon / (1 - a**2 / b**2)
    q_s_b_expr = -4 * pi * L * V * epsilon * a**2 / (b**2 * (1 - a**2 / b**2))

    # Print the results
    print("Total volume charge (q_v):")
    sympy.pprint(q_v_expr)
    print("\nTotal surface charge on inner electrode (q_s(r=a)):")
    sympy.pprint(q_s_a_expr)
    print("\nTotal surface charge on outer electrode (q_s(r=b)):")
    sympy.pprint(q_s_b_expr)

solve_charges()