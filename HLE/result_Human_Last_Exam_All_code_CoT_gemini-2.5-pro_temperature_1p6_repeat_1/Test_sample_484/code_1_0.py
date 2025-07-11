import sympy

def solve_synaptic_plasticity():
    """
    Performs a steady-state analysis of the given biophysical model
    to derive an expression for the change in synaptic efficacy.
    """
    # Define the symbolic variables from the model
    # w_i is the synaptic efficacy (W_i)
    # v_i is the presynaptic accumulator (M_i)
    # u_i is the postsynaptic accumulator (Y)
    w_i, v_i, u_i = sympy.symbols('w_i v_i u_i')
    
    # Define the constants
    alpha, beta, eta, tau_W = sympy.symbols('alpha beta eta tau_W')
    
    # Define placeholders for the fast variables P_i and B_i
    P_i, B_i = sympy.symbols('P_i B_i')

    # Step 1: Write down the steady-state equations for P_i and B_i.
    # Original equations:
    # tau_P * dP_i/dt = -P_i + (1-eta)*Y - M_i*P_i
    # tau_P * dB_i/dt = -B_i + eta*Y + M_i*P_i
    # At steady-state (d/dt = 0) and with substitutions (Y->u_i, M_i->v_i):
    # 0 = -P_i + (1-eta)*u_i - v_i*P_i  => P_i(1 + v_i) = (1-eta)*u_i
    # 0 = -B_i + eta*u_i + v_i*P_i      => B_i = eta*u_i + v_i*P_i
    
    # Solve for P_i at steady state
    P_i_ss = sympy.solve(sympy.Eq(-P_i + (1 - eta) * u_i - v_i * P_i, 0), P_i)[0]
    
    # Solve for B_i at steady state by substituting P_i_ss
    B_i_ss = eta * u_i + v_i * P_i_ss
    
    # Step 2: Substitute the steady-state expressions into the equation for w_i dynamics.
    # Original equation: tau_W * dW_i/dt = alpha*P_i + beta*B_i
    # We want an expression for the RHS.
    
    tau_W_dwidt = alpha * P_i_ss + beta * B_i_ss
    
    # Step 3: Simplify the expression.
    tau_W_dwidt_simplified = sympy.simplify(tau_W_dwidt)
    
    # Let's collect terms with respect to u_i for better readability
    tau_W_dwidt_collected = sympy.collect(tau_W_dwidt_simplified, u_i)
    
    # Step 4: Define the constant rho and substitute it into the expression.
    # The term in the parenthesis is (alpha*(1-eta) + beta*eta + beta*v_i).
    # We can define rho such that the expression simplifies. Let's define the plasticity
    # threshold rho as the value of v_i for which the net change is zero.
    # alpha*(1-eta) + beta*eta + beta*v_i = 0 => v_i = - (alpha*(1-eta) + beta*eta) / beta
    
    rho = sympy.Symbol('rho')
    # From the definition rho = -(alpha*(1-eta) + beta*eta) / beta, we have:
    # alpha*(1-eta) + beta*eta = -beta*rho
    
    # The numerator of our expression is u_i * (alpha - alpha*eta + beta*eta + beta*v_i)
    # This is u_i * ( (alpha*(1-eta) + beta*eta) + beta*v_i )
    # Substitute the definition of rho: u_i * ( -beta*rho + beta*v_i )
    # This simplifies to u_i * beta * (v_i - rho)
    
    final_numerator = u_i * beta * (v_i - rho)
    final_denominator = 1 + v_i
    final_expression = final_numerator / final_denominator

    # Final result
    print("The derived expression for the rate of change of synaptic efficacy is:")
    lhs = sympy.Symbol('tau_W*dw_i/dt')
    final_equation = sympy.Eq(lhs, final_expression)
    print(final_equation)
    
    # Define rho explicitly
    rho_definition = sympy.Eq(rho, - (alpha * (1 - eta) + beta * eta) / beta)
    print("\nWhere the constant rho is defined as:")
    print(rho_definition)

solve_synaptic_plasticity()