import sympy as sp

def derive_simplified_plasticity_rule():
    """
    Performs a steady-state analysis of the given biophysical model to derive
    a simplified synaptic plasticity rule using SymPy.
    """
    print("Step 1: Define all symbolic variables from the model.")
    # Define time constants and model parameters
    tau_W, phi, eta, alpha, beta = sp.symbols('tau_W phi eta alpha beta', real=True, positive=True)
    # Define dynamic variables. x_i is the presynaptic input rate.
    Y, M_i, P_i, B_i, x_i = sp.symbols('Y M_i P_i B_i x_i')
    # Y_sum represents the weighted sum of all inputs: sum_j(w_j*x_j)
    Y_sum = sp.Symbol('Y_sum')
    # w_dot_i represents the time derivative of the synaptic weight w_i
    w_dot_i = sp.Symbol('w_dot_i')
    print("-" * 30)

    print("Step 2: Solve for the steady-state of fast variables (M, Y, P, B).")
    # Steady state means the time derivatives are zero.

    # For M_i from Eq1: 0 = -M_i + phi * x_i
    ss_M_i = phi * x_i
    print("\nSteady-state of presynaptic MMP9 (M_i):")
    sp.pprint(sp.Eq(M_i, ss_M_i))

    # For Y: 0 = -Y + sum(w_j*x_j)
    ss_Y = Y_sum
    print("\nSteady-state of postsynaptic Calcium (Y):")
    sp.pprint(sp.Eq(Y, ss_Y))

    # For P_i and B_i, we solve the system:
    # 0 = -P_i + (1-eta)*Y - M_i*P_i
    # 0 = -B_i + eta*Y + M_i*P_i
    # A key insight from summing these equations is that at steady state, P_i + B_i = Y.
    # From the first equation: P_i * (1 + M_i) = (1 - eta) * Y
    ss_P_i = (1 - eta) * Y / (1 + M_i)
    ss_B_i = Y - ss_P_i  # Using the P_i + B_i = Y relationship
    print("\nSteady-state of proBDNF (P_i):")
    sp.pprint(sp.Eq(P_i, ss_P_i))
    print("\nSteady-state of BDNF (B_i):")
    sp.pprint(sp.Eq(B_i, ss_B_i))
    print("-" * 30)
    
    print("Step 3: Substitute steady-state expressions into the weight dynamics equation.")
    # Original equation: tau_W * dw_i/dt = alpha*P_i + beta*B_i
    w_dot_i_expr = alpha * ss_P_i + beta * ss_B_i

    # Now substitute the expressions for M_i and Y
    w_dot_i_final_expr = w_dot_i_expr.subs({
        Y: ss_Y,
        M_i: ss_M_i
    })
    
    # Simplify the expression
    w_dot_i_final_expr = sp.simplify(w_dot_i_final_expr)
    print("\nExpression for weight dynamics (tau_W * dw_i/dt):")
    sp.pprint(w_dot_i_final_expr)
    print("-" * 30)

    print("Step 4: Introduce simplified variables u, v_i, and rho.")
    u, v_i, rho = sp.symbols('u v_i rho')

    # Define the new variables:
    # Presynaptic accumulator v_i is the steady-state MMP9 level.
    # Postsynaptic accumulator u is the steady-state calcium level.
    # rho is the ratio of LTD/LTP strengths.
    print("\nDefinitions of new variables:")
    print(f"  v_i = M_i = {sp.sstr(ss_M_i, full_prec=False)}")
    print(f"  u   = Y   = Y_sum (i.e., sum over j of w_j*x_j)")
    print(f"  rho = alpha / beta")
    
    # Perform the substitutions to get the final form
    # Replace Y_sum with u, and phi*x_i with v_i
    final_eq_rhs = w_dot_i_final_expr.subs({
        Y_sum: u,
        phi * x_i: v_i
    })

    # Replace alpha with rho*beta
    final_eq_rhs = final_eq_rhs.subs({alpha: rho * beta})

    # Factor and simplify to the desired form: u*beta*(1 + ...)
    final_eq_rhs = sp.factor(final_eq_rhs, u, beta)
    term_in_paren = sp.simplify(final_eq_rhs / (u * beta))
    
    # Use partial fraction decomposition to get the "1 + ..." form
    term_in_paren = sp.apart(term_in_paren, v_i)

    final_eq_rhs_nice = u * beta * term_in_paren
    final_eq_str = f"tau_W * dw_i/dt = {sp.sstr(final_eq_rhs_nice, full_prec=False)}"
    
    print("\nFinal derived expression:")
    print("The final equation for the change in synaptic efficacy is:")
    # We print each component part of the equation's right hand side
    print(f"  Postsynaptic accumulator part: {u}")
    print(f"  LTP strength part: {beta}")
    # The term_in_paren is the complex multiplier
    multiplier = sp.sstr(term_in_paren, full_prec=False)
    print(f"  Modulatory part: {multiplier}")
    
    final_equation_obj = sp.Eq(tau_W * w_dot_i, final_eq_rhs_nice)
    print("\nFull equation:")
    sp.pprint(final_equation_obj)

if __name__ == '__main__':
    derive_simplified_plasticity_rule()