import sympy

def display_steady_state_equation():
    """
    This function presents the derived steady-state equation for synaptic plasticity
    from the biophysical model.

    It uses the sympy library to construct and display the final mathematical
    expressions for the synaptic efficacy dynamics and the definition of the
    constant rho.
    """
    # Define the symbolic variables required for the final expressions.
    # w_i_dot represents the time derivative of w_i.
    tau_w, w_i_dot = sympy.symbols('tau_w \\dot{w_i}')
    u_i, v_i = sympy.symbols('u_i v_i')
    beta, rho, alpha, eta = sympy.symbols('beta, rho, alpha, eta')

    # Based on the steady-state derivation, we construct the final equation.
    # The left-hand side (LHS) of the equation is tau_w * d(w_i)/dt.
    lhs = tau_w * w_i_dot

    # The right-hand side (RHS) is the derived expression u_i * (beta - rho / (1 + v_i)).
    rhs = u_i * (beta - rho / (1 + v_i))

    # We create a symbolic equation object for pretty printing.
    final_equation = sympy.Eq(lhs, rhs)

    # We also define the expression for the constant rho, which groups
    # the parameters alpha, beta, and eta.
    rho_definition = sympy.Eq(rho, (beta - alpha) * (1 - eta))

    # Print the final results in a clear format.
    print("The derived steady-state expression for synaptic efficacy is:")
    sympy.pprint(final_equation, use_unicode=True)

    print("\nwhere the variables are defined as:")
    print(" w_i: Synaptic efficacy")
    print(" u_i: Postsynaptic accumulator (total weighted input, y)")
    print(" v_i: Presynaptic accumulator (presynaptic MMP9 level, m_i)")

    print("\nAnd the constant rho is defined as:")
    sympy.pprint(rho_definition, use_unicode=True)


if __name__ == '__main__':
    display_steady_state_equation()