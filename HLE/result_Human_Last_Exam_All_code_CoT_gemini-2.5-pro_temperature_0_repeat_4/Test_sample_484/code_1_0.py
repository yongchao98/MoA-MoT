import sympy

def derive_synaptic_plasticity_rule():
    """
    This function performs a steady-state analysis of the given biophysical model
    to derive a simplified synaptic plasticity rule.
    It prints the derivation steps and the final result.
    """

    # Define symbols
    # Using sympy for clarity, though not strictly necessary for this derivation
    tau_w, w_i, t = sympy.symbols('tau_w w_i t')
    alpha, beta, eta, phi = sympy.symbols('alpha beta eta phi')
    M_i, Y, P_i, B_i = sympy.symbols('M_i Y P_i B_i')
    x_i = sympy.Symbol('x_i')
    w_j, x_j = sympy.symbols('w_j x_j')
    N = sympy.Symbol('N')
    u, v_i, rho = sympy.symbols('u v_i rho')

    print("Derivation of the simplified synaptic plasticity rule:")
    print("="*50)

    # Step 1: Steady-state analysis for M_i and Y
    print("Step 1: Find the steady-state values for the fast variables M_i and Y.")
    print("The dynamics are assumed to be much faster than w_i, so we set their time derivatives to 0.")
    
    # M_i steady state
    eq1 = -M_i + phi * x_i
    M_i_ss = sympy.solve(eq1, M_i)[0]
    print(f"\nFrom tau_M * dM_i/dt = -M_i + phi*x_i = 0, we get:")
    print(f"M_i_ss = {M_i_ss}")

    # Y steady state
    # The term w^T(t)X(t) is the sum over all synapses j=1 to N of w_j * x_j
    sum_w_x = sympy.Symbol('sum_j(w_j*x_j)')
    eq2 = -Y + sum_w_x
    Y_ss = sympy.solve(eq2, Y)[0]
    print(f"\nFrom tau_Y * dY/dt = -Y + sum_j(w_j*x_j) = 0, we get:")
    print(f"Y_ss = {Y_ss}")
    print("-" * 50)

    # Step 2: Define u and v_i
    print("Step 2: Define the postsynaptic accumulator 'u' and presynaptic accumulator 'v_i'.")
    # u is the postsynaptic accumulator, representing steady-state calcium
    u_def = Y_ss
    print(f"\nWe define the postsynaptic accumulator u as the steady-state value of Y:")
    print(f"u = {u_def}")
    # v_i is the presynaptic accumulator, representing steady-state MMP9
    v_i_def = M_i_ss
    print(f"\nWe define the presynaptic accumulator v_i as the steady-state value of M_i:")
    print(f"v_i = {v_i_def}")
    print("-" * 50)

    # Step 3: Steady-state analysis for P_i and B_i
    print("Step 3: Find the steady-state values for P_i and B_i in terms of u and v_i.")
    
    # P_i steady state
    eq3 = -P_i + (1 - eta) * u - v_i * P_i
    P_i_ss = sympy.solve(eq3, P_i)[0]
    print(f"\nFrom tau_P * dP_i/dt = -P_i + (1-eta)*Y - M_i*P_i = 0, we substitute Y=u and M_i=v_i:")
    print(f"-P_i + (1-eta)*u - v_i*P_i = 0  =>  P_i*(1 + v_i) = (1-eta)*u")
    print(f"P_i_ss = {P_i_ss}")

    # B_i steady state
    eq4 = -B_i + eta * u + v_i * P_i_ss
    B_i_ss = sympy.solve(eq4, B_i)[0]
    print(f"\nFrom tau_P * dB_i/dt = -B_i + eta*Y + M_i*P_i = 0, we substitute Y=u, M_i=v_i, and P_i_ss:")
    print(f"B_i_ss = eta*u + v_i * ({P_i_ss})")
    B_i_ss_simplified = sympy.simplify(B_i_ss)
    print(f"B_i_ss = {B_i_ss_simplified}")
    print("-" * 50)

    # Step 4: Substitute into the equation for w_i
    print("Step 4: Substitute the steady-state expressions into the equation for w_i.")
    dw_dt_expr = alpha * P_i_ss + beta * B_i_ss
    print(f"\nThe equation for w_i is: tau_w * dw_i/dt = alpha*P_i + beta*B_i")
    print(f"Substituting P_i_ss and B_i_ss:")
    print(f"tau_w * dw_i/dt = alpha*({P_i_ss}) + beta*({B_i_ss_simplified})")
    
    dw_dt_simplified = sympy.simplify(dw_dt_expr)
    print(f"\nCombining terms, we get:")
    print(f"tau_w * dw_i/dt = {dw_dt_simplified}")
    print("-" * 50)

    # Step 5: Final simplification and definition of rho
    print("Step 5: Rearrange the expression and define the constant 'rho'.")
    # Factor out beta*u / (1+v_i)
    term_in_brackets = sympy.expand(dw_dt_simplified * (1 + v_i) / (beta * u))
    print(f"\nFactoring out (beta*u / (1+v_i)), we get:")
    print(f"tau_w * dw_i/dt = (beta*u / (1+v_i)) * ({term_in_brackets})")
    
    # Define rho
    rho_def_expr = -(term_in_brackets - v_i)
    print(f"\nWe define a constant rho such that the term in brackets is (v_i - rho).")
    print(f"rho = -({term_in_brackets - v_i})")
    rho_def_simplified = sympy.simplify(rho_def_expr)
    print(f"rho = {rho_def_simplified}")

    final_expression = f"tau_w*dw_i/dt = (beta * u / (1 + v_i)) * (v_i - rho)"
    
    print("\n" + "="*50)
    print("Final Result:")
    print("The derived expression for the change in synaptic efficacy is:")
    print(f"\n{final_expression}\n")
    print("where:")
    print(f"  u = {u_def} (postsynaptic activity)")
    print(f"  v_i = {v_i_def} (presynaptic activity)")
    print(f"  rho = {rho_def_simplified} (plasticity threshold)")
    print("="*50)
    
    # Final answer in the required format
    final_answer_string = f"tau_w*dw_i/dt = (beta * u / (1 + v_i)) * (v_i - rho)"
    print(f"\n<<<{final_answer_string}>>>")

if __name__ == '__main__':
    derive_synaptic_plasticity_rule()