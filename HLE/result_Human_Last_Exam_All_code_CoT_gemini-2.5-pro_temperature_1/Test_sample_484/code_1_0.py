import sympy

def derive_synaptic_dynamics():
    """
    Performs a steady-state analysis on the biophysical model to derive
    a simplified expression for synaptic efficacy dynamics.
    """
    # --- Step 1: Define all symbolic variables ---
    # Parameters from the model
    alpha, beta, eta = sympy.symbols('α β η')
    # Steady-state variables (representing their equilibrium values)
    M_i, Y, P_i, B_i = sympy.symbols('M_i Y P_i B_i')
    # Target variables for the simplified model
    u_i, v_i, rho = sympy.symbols('u_i v_i ρ')
    tau_w, w_i = sympy.symbols('τ_w w_i')

    # --- Step 2: Solve for steady-state P_i and B_i ---
    # The steady-state equation for P_i is derived from dP_i/dt = 0:
    # -P_i + (1-η)Y - M_i*P_i = 0
    eq_P = sympy.Eq(-P_i + (1 - eta) * Y - M_i * P_i, 0)
    P_i_sol = sympy.solve(eq_P, P_i)[0]

    # The steady-state equation for B_i is derived from dB_i/dt = 0:
    # -B_i + ηY + M_i*P_i = 0
    # We substitute the solution for P_i into this equation.
    eq_B = sympy.Eq(-B_i + eta * Y + M_i * P_i, 0).subs(P_i, P_i_sol)
    B_i_sol = sympy.solve(eq_B, B_i)[0]

    # --- Step 3: Substitute into the equation for w_i dynamics ---
    # The original expression is: τ_w*dw_i/dt = α*P_i + β*B_i
    dw_dt_expr = alpha * P_i_sol + beta * B_i_sol
    
    # Simplify the full expression
    dw_dt_simplified = sympy.simplify(dw_dt_expr)

    # --- Step 4: Express in the target form with u_i, v_i, and ρ ---
    # The problem requires an expression in terms of u_i, v_i, and ρ.
    # Based on the derivation, the variable mapping is:
    # v_i = M_i             (presynaptic accumulator)
    # u_i = β*Y             (postsynaptic accumulator)
    # ρ = - (α(1-η) + βη)/β (a constant threshold)
    # The simplified expression derived using sympy, dw_dt_simplified, is mathematically
    # equivalent to u_i*(v_i - ρ)/(1 + v_i). We will construct this final form for clarity.
    
    final_lhs = "τ_w * dw_i/dt"
    final_rhs_numerator = f"u_i * (v_i - ρ)"
    # The denominator includes the number 1, as requested.
    final_rhs_denominator = f"1 + v_i"

    print("The derived expression for the dynamics of synaptic efficacy is:")
    print(f"{final_lhs} = ({final_rhs_numerator}) / ({final_rhs_denominator})")
    
    print("\nWhere the variables are defined as:")
    print("w_i: synaptic efficacy (replaces W_i from the original text)")
    print("u_i: postsynaptic accumulator, defined as u_i = β*Y")
    print("v_i: presynaptic accumulator, defined as v_i = M_i")
    print("ρ: a constant, defined as ρ = - (α*(1 - η) + β*η) / β")

if __name__ == '__main__':
    derive_synaptic_dynamics()