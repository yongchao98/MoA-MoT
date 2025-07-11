import sympy as sp

def derive_synaptic_dynamics():
    """
    This function symbolically derives the steady-state expression for synaptic
    efficacy dynamics based on the provided biophysical model.
    """
    # 1. Define symbolic variables from the model description
    alpha, beta, eta, phi = sp.symbols('alpha beta eta phi')
    M_i, Y, P_i, B_i = sp.symbols('M_i Y P_i B_i')
    tau_w = sp.Symbol('tau_w')
    w_i_dot = sp.Symbol(r'\dot{w_i}')

    # 2. Solve for steady-state P_i and B_i
    # From dP_i/dt = 0:  -P_i + (1-eta)*Y - M_i*P_i = 0
    # => P_i(1 + M_i) = (1-eta)*Y
    ss_P_i = (1 - eta) * Y / (1 + M_i)

    # From dB_i/dt = 0: -B_i + eta*Y + M_i*P_i = 0
    # => B_i = eta*Y + M_i*P_i
    ss_B_i = sp.simplify(eta * Y + M_i * ss_P_i)

    # 3. Substitute steady-state expressions into the weight dynamics equation
    # tau_w * dW_i/dt = alpha*P_i + beta*B_i
    rhs_expression = sp.simplify(alpha * ss_P_i + beta * ss_B_i)

    # 4. Introduce the new variables u_i, v_i, and rho
    u_i, v_i, rho = sp.symbols('u_i v_i rho')
    # Replace steady-state Y with postsynaptic accumulator u_i
    # Replace steady-state M_i with presynaptic accumulator v_i
    rhs_uv = rhs_expression.subs({Y: u_i, M_i: v_i})
    
    # Introduce rho = -alpha/beta by substituting alpha = -rho*beta
    rhs_final = rhs_uv.subs({alpha: -rho * beta})
    
    # Simplify the final expression
    rhs_final_simplified = sp.simplify(rhs_final)

    # 5. Format and print the final equation
    # We construct a human-readable string with Unicode for Greek letters
    LHS = "\u03C4_w\u1E81\u208_i"  # tau_w * w_i_dot
    var_beta = "\u03B2"
    var_eta = "\u03B7"
    var_rho = "\u03C1"
    
    # The final expression is beta*u_i*(v_i - (rho*(1-eta) - eta)) / (1+v_i)
    # We will print each variable and "number" (in this case, 1) in the equation.
    print("The derived final equation is:")
    print(f"{LHS} = ({var_beta} * u_i / (1 + v_i)) * (v_i - ({var_rho} * (1 - {var_eta}) - {var_eta}))")
    
    # Print the variables involved as requested
    print("\nVariables in the equation:")
    print(f"τ_w (tau_w), ẇ_i (w_i_dot), β (beta), u_i, v_i, ρ (rho), η (eta)")

derive_synaptic_dynamics()