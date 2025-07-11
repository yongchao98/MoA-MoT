import sympy
from sympy import symbols, pi, pprint

def solve_electrostatics_problem():
    """
    This function prints the symbolic expressions for the total volume charge and
    the surface charges on the electrodes based on the selected answer choice.
    """
    # Define symbolic variables
    V, epsilon, pi, L, a, b = symbols('V varepsilon pi L a b')

    # Formulas from Option B
    q_v_expr = -4 * V * epsilon * pi * L / (1 - a**2/b**2)
    q_s_a_expr = 2 * pi * L * V * epsilon / (1 - a**2/b**2)
    q_s_b_expr = -4 * pi * L * V * epsilon * a**2 / (b**2 * (1 - a**2/b**2))

    print("The expressions for the charges are:")
    print("\nTotal volume charge (q_v):")
    pprint(q_v_expr, use_unicode=True)
    
    print("\nTotal surface charge on inner electrode q_s(r=a):")
    pprint(q_s_a_expr, use_unicode=True)
    
    print("\nTotal surface charge on outer electrode q_s(r=b):")
    pprint(q_s_b_expr, use_unicode=True)

solve_electrostatics_problem()