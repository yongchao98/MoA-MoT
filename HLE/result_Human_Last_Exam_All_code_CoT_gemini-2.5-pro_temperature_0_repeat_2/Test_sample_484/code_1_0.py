import sympy

def derive_simplified_model():
    """
    Performs a steady-state analysis of the biophysical model using sympy
    to derive the simplified expression for synaptic efficacy dynamics.
    """
    # Define the parameters and variables of the model as symbolic objects
    alpha, beta, eta, phi = sympy.symbols('alpha beta eta phi')
    tau_w = sympy.symbols('tau_W')
    nu_i, w_i = sympy.symbols('nu_i w_i')
    
    # Steady-state variables (represented with _ss)
    M_i_ss, Y_ss = sympy.symbols('M_i_ss Y_ss')
    
    # Target variables for the simplified system
    u_i, v_i, rho = sympy.symbols('u_i v_i rho')
    
    # The derivative of w_i, denoted as w_i_dot
    w_i_dot = sympy.Symbol(r'\dot{w_i}')

    # --- Step 1: Express steady-state variables ---
    # From dM/dt = 0, M_i = phi * nu_i. This is our presynaptic accumulator, v_i.
    M_i_ss_expr = phi * nu_i
    
    # From dY/dt = 0, Y = sum(w_j * nu_j). This is our postsynaptic accumulator, u_i.
    # We represent it symbolically as u_i as requested.
    Y_ss_expr = u_i
    
    # From dP/dt = 0 and dB/dt = 0, we solve for P_i and B_i
    # P_i_ss * (1 + M_i_ss) = (1 - eta) * Y_ss
    P_i_ss_expr = (1 - eta) * Y_ss_expr / (1 + M_i_ss_expr)
    
    # B_i_ss = eta * Y_ss + M_i_ss * P_i_ss
    B_i_ss_expr = eta * Y_ss_expr + M_i_ss_expr * P_i_ss_expr
    B_i_ss_expr = sympy.simplify(B_i_ss_expr)

    # --- Step 2: Substitute into the weight dynamics equation ---
    # tau_W * dW_i/dt = alpha * P_i + beta * B_i
    # We use the steady-state expressions for P_i and B_i
    rhs_expr = alpha * P_i_ss_expr + beta * B_i_ss_expr
    rhs_expr = sympy.simplify(rhs_expr)

    # --- Step 3: Introduce the new variables and simplify ---
    # The expression is currently in terms of u_i and (phi*nu_i).
    # We define v_i = phi * nu_i (presynaptic accumulator)
    rhs_expr = rhs_expr.subs(M_i_ss_expr, v_i)
    
    # The expression is: u_i*(alpha*(1 - eta) + beta*(eta + v_i))/(v_i + 1)
    # We can rearrange this to: u_i * ( (alpha*(1-eta) + beta*eta) + beta*v_i ) / (v_i + 1)
    # Let's define rho such that beta*rho = alpha*(1-eta) + beta*eta
    rho_definition = (alpha*(1-eta) + beta*eta) / beta
    
    # Now we create the final expression using rho
    final_rhs = beta * u_i * (rho + v_i) / (1 + v_i)

    # --- Step 4: Print the final results ---
    print("The steady-state analysis reduces the system to the following expression for the change in synaptic efficacy:")
    
    # Create the full equation for printing
    lhs = tau_w * w_i_dot
    final_equation = sympy.Eq(lhs, final_rhs)
    
    print("\nFinal Equation:")
    sympy.pprint(final_equation, use_unicode=True)
    
    print("\nWhere the variables are defined as:")
    print(f"  v_i = M_i_ss = phi * nu_i  (presynaptic accumulator)")
    print(f"  u_i = Y_ss = sum over j of (w_j * nu_j)  (postsynaptic accumulator)")
    
    print("\nAnd the constant rho is defined as:")
    sympy.pprint(sympy.Eq(rho, rho_definition), use_unicode=True)

    print("\nAs requested, here is each number in the final equation:")
    # The only explicit number in the final simplified form is '1' in the denominator.
    print("1")

if __name__ == '__main__':
    derive_simplified_model()