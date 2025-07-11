import sympy as sp

def print_solution():
    """
    This function prints the formulas for total volume charge and surface charges as given in the selected answer choice.
    """
    # Define symbolic variables for clarity in the formulas
    V, epsilon, pi, L, a, b = sp.symbols('V varepsilon pi L a b')

    # Formulas from Option B
    q_v_expr_str = "Total volume charge = q_v = -4 * V * varepsilon * pi * L / (1 - a**2/b**2)"
    q_s_a_expr_str = "Total surface charge on inner electrode = q_s(r = a) = 2 * pi * L * V * varepsilon / (1 - a**2/b**2)"
    q_s_b_expr_str = "Total surface charge on outer electrode = q_s(r = b) = -4 * pi * L * V * varepsilon * a**2 / (b**2 * (1 - a**2/b**2))"

    print(q_v_expr_str)
    print(q_s_a_expr_str)
    print(q_s_b_expr_str)

print_solution()