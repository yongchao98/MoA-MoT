import sympy

def derive_and_print_equation():
    """
    This function derives the steady-state equation for synaptic efficacy and prints the results.
    The derivation steps are outlined in the text above.
    """

    # Define symbols for clarity in the output.
    # We use sympy for pretty printing, but the derivation is done manually.
    tau_W, w_i_dot = sympy.symbols('tau_W, dot(w_i)')
    alpha, beta, eta, phi = sympy.symbols('alpha, beta, eta, phi')
    u_i, v_i, rho = sympy.symbols('u_i, v_i, rho')
    w_j, x_j, x_i = sympy.symbols('w_j, x_j, x_i')
    j = sympy.Idx('j')

    # Define the accumulators and the constant rho based on the derivation.
    # The problem asks for a smaller system with u_i and v_i.
    # Based on the steady-state analysis:
    # v_i represents the presynaptic state (M_i).
    # u_i represents the postsynaptic state (Y).
    # rho combines the plasticity parameters.

    v_i_def = phi * x_i
    u_i_def = sympy.Sum(w_j * x_j, j)
    rho_def = eta + (alpha / beta) * (1 - eta)

    # Construct the final equation for the dynamics of the synaptic efficacy w_i.
    # The derived expression is: tau_W * dot(w_i) = beta * u_i * (v_i + rho) / (1 + v_i)

    final_equation_lhs = tau_W * w_i_dot
    final_equation_rhs = beta * u_i * (v_i + rho) / (1 + v_i)
    final_equation = sympy.Eq(final_equation_lhs, final_equation_rhs)

    # Print the definitions and the final equation.
    print("This script presents the result of the steady-state analysis.")
    print("------------------------------------------------------------\n")

    print("Definition of the presynaptic accumulator (v_i):")
    sympy.pprint(sympy.Eq(v_i, v_i_def), use_unicode=False)
    print("\nThis represents the steady-state level of presynaptic MMP9, proportional to the presynaptic input x_i.\n")

    print("Definition of the postsynaptic accumulator (u_i):")
    print("# Note: The postsynaptic accumulator Y is common to all synapses, so u_i is independent of i.")
    sympy.pprint(sympy.Eq(u_i, u_i_def), use_unicode=False)
    print("\nThis represents the steady-state postsynaptic calcium level, proportional to the weighted sum of all inputs.\n")

    print("Definition of the constant rho:")
    sympy.pprint(sympy.Eq(rho, rho_def), use_unicode=False)
    print("\nThis constant combines the parameters for LTD (alpha), LTP (beta), and the baseline BDNF/proBDNF ratio (eta).\n")

    print("------------------------------------------------------------")
    print("Final simplified expression for the synaptic efficacy dynamics:")
    print("------------------------------------------------------------")
    sympy.pprint(final_equation, use_unicode=False)

    # To fulfill the "output each number in the final equation" request,
    # we can break down and print the components of the final equation.
    print("\n--- Components of the Final Equation ---")
    print(f"1. Rate of change term: {sympy.pretty(final_equation_lhs, use_unicode=False)}")
    print(f"2. Postsynaptic drive term: {sympy.pretty(beta * u_i, use_unicode=False)}")
    print(f"3. Presynaptic modulation term: {sympy.pretty((v_i + rho) / (1 + v_i), use_unicode=False)}")


if __name__ == '__main__':
    derive_and_print_equation()
    print("\n<<<tau_W*dot(w_i) = beta*u_i*(v_i + rho)/(1 + v_i)>>>")
