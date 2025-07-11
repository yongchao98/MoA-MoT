import sympy

def derive_simplified_model():
    """
    Performs a steady-state analysis on the biophysical model using sympy
    to derive the simplified learning rule.
    """
    # 1. Define all symbols for the derivation.
    # We use names without subscripts for easier handling in sympy.
    # Base parameters
    alpha, beta, eta, phi = sympy.symbols('alpha beta eta phi')

    # State variables (steady-state) and inputs
    u, v = sympy.symbols('u_i v_i') # Post- and Pre-synaptic accumulators
    w, nu = sympy.symbols('w_j nu_j') # Generic weight and firing rate for the sum
    
    # Target expression symbols
    tau_W, w_i_dot = sympy.symbols('tau_W dot_w_i')
    rho = sympy.Symbol('rho')

    print("--- Step 1: Defining variables based on steady-state analysis ---")
    
    # Presynaptic accumulator v_i is the steady-state of M_i
    # M_i_ss = phi * nu_i. We define v_i = M_i_ss.
    print("Presynaptic accumulator (v_i) is the steady-state MMP9 level: v_i = phi * nu_i")
    
    # Postsynaptic accumulator u_i is the steady-state of Y
    # Y_ss = sum(w_j * nu_j). We define u_i = Y_ss.
    print("Postsynaptic accumulator (u_i) is the steady-state Calcium level: u_i = Sum(w_j * nu_j)")
    
    print("\n--- Step 2: Expressing steady-state P_i and B_i in terms of u_i and v_i ---")

    # From dP_i/dt = 0: P_i_ss = (1-eta)*Y_ss / (1+M_i_ss)
    Pi_ss = (1 - eta) * u / (1 + v)
    
    # From dB_i/dt = 0: B_i_ss = (eta*Y_ss + M_i_ss*P_i_ss)
    # This simplifies to B_i_ss = Y_ss * (eta + M_i_ss) / (1 + M_i_ss)
    Bi_ss = u * (eta + v) / (1 + v)
    
    print("Steady-state proBDNF (P_i_ss):")
    sympy.pprint(sympy.Eq(sympy.Symbol("P_i_ss"), Pi_ss))
    print("\nSteady-state BDNF (B_i_ss):")
    sympy.pprint(sympy.Eq(sympy.Symbol("B_i_ss"), Bi_ss))

    print("\n--- Step 3: Substituting into the synaptic efficacy equation ---")
    # tau_W * dW_i/dt = alpha*P_i + beta*B_i
    learning_rule_expr = alpha * Pi_ss + beta * Bi_ss
    
    print("Initial expression for tau_W * dw_i/dt:")
    sympy.pprint(learning_rule_expr)

    print("\n--- Step 4: Simplifying and refactoring to the target form ---")
    # Simplify the expression by combining terms over a common denominator
    simplified_expr = sympy.simplify(learning_rule_expr)
    
    # The simplified expression is u_i * (alpha*(1-eta) + beta*eta + beta*v_i) / (1+v_i)
    # We can refactor this into the form u_i * beta * (v_i - rho) / (1+v_i)
    
    # Let's define rho such that this refactoring works.
    # We want: alpha*(1-eta) + beta*eta + beta*v_i = beta * (v_i - rho)
    # alpha - alpha*eta + beta*eta + beta*v_i = beta*v_i - beta*rho
    # alpha - alpha*eta + beta*eta = -beta*rho
    # rho = -(alpha - alpha*eta + beta*eta) / beta
    rho_definition = -(alpha * (1 - eta) + beta * eta) / beta
    
    # Create the final target expression using the symbol 'rho'
    final_rhs = u * beta * (v - rho) / (1 + v)
    final_equation = sympy.Eq(tau_W * w_i_dot, final_rhs)
    
    print("\nThe final derived expression for the learning rule is:")
    sympy.pprint(final_equation, use_unicode=True)
    
    print("\nWhere the composite constant rho is defined as:")
    sympy.pprint(sympy.Eq(rho, rho_definition), use_unicode=True)

    # Final check: Ensure the derived expression is identical to the refactored one
    check_expr = simplified_expr - final_rhs.subs(rho, rho_definition)
    assert sympy.simplify(check_expr) == 0

    print("\n--- Final Equation Components ---")
    print("The numbers and symbols (parameters and variables) in the final equation are:")
    # The .free_symbols attribute gives all the individual symbols in an expression
    all_symbols = final_equation.free_symbols.union(rho_definition.free_symbols)
    # Sort for consistent output
    sorted_symbols = sorted(list(all_symbols), key=lambda s: s.name)
    for s in sorted_symbols:
        print(s.name)


if __name__ == '__main__':
    derive_simplified_model()