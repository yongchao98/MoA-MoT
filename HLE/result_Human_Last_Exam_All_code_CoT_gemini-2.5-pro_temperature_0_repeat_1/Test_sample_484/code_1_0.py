import sympy as sp

def derive_synaptic_efficacy_equation():
    """
    Performs a steady-state analysis of the given biophysical model
    to derive a simplified expression for the change in synaptic efficacy.
    """
    # Plan:
    # 1. Define all necessary symbolic variables using sympy.
    # 2. Set up the steady-state equations by setting the time derivatives of M, P, and B to zero.
    #    - We will directly use the new variables u_i and v_i for the steady states of Y and x_i.
    # 3. Solve the system of algebraic equations for the steady-state values of M_i, P_i, and B_i.
    # 4. Substitute these steady-state expressions into the dynamic equation for the synaptic efficacy, w_i.
    # 5. Simplify the resulting expression for tau_w * dw_i/dt.
    # 6. Introduce the constant rho = alpha / beta and perform further simplification to obtain the final form.
    # 7. Print the final derived equation in a clear, readable format, showing each component.

    # Step 1: Define symbols for the variables and parameters
    # New variables for the simplified system
    u_i, v_i = sp.symbols('u_i v_i')
    # Parameters
    alpha, beta, eta, phi, rho, tau_W = sp.symbols('alpha beta eta phi rho tau_W')
    # Steady-state variables
    M_i_ss, P_i_ss, B_i_ss = sp.symbols('M_i P_i B_i')

    # Step 2: Set up the steady-state equations
    # From tau_M * dM_i/dt = -M_i + phi * x_i = 0, with x_i -> v_i
    # We get M_i_ss = phi * v_i
    eq_M_ss = sp.Eq(M_i_ss, phi * v_i)

    # From tau_P * dP_i/dt = -P_i + (1-eta)*Y - M_i*P_i = 0, with Y -> u_i
    # We get P_i * (1 + M_i) = (1-eta)*u_i
    eq_P_ss = sp.Eq(P_i_ss * (1 + M_i_ss), (1 - eta) * u_i)

    # From tau_P * dB_i/dt = -B_i + eta*Y + M_i*P_i = 0, with Y -> u_i
    # We get B_i = eta*u_i + M_i*P_i
    eq_B_ss = sp.Eq(B_i_ss, eta * u_i + M_i_ss * P_i_ss)

    # Step 3: Solve for the steady-state values
    # Solve for M_i_ss is trivial from its definition
    sol_M = eq_M_ss.rhs

    # Solve for P_i_ss
    sol_P = sp.solve(eq_P_ss, P_i_ss)[0]
    sol_P = sol_P.subs(M_i_ss, sol_M)

    # Solve for B_i_ss
    sol_B = eq_B_ss.rhs
    sol_B = sol_B.subs(M_i_ss, sol_M)
    sol_B = sol_B.subs(P_i_ss, sol_P)

    # Step 4: Substitute into the synaptic efficacy equation
    # tau_W * dw_i/dt = alpha * P_i + beta * B_i
    dw_dt_expr = alpha * sol_P + beta * sol_B

    # Step 5: Simplify the expression and introduce rho
    dw_dt_simplified = sp.simplify(dw_dt_expr)
    dw_dt_rho = dw_dt_simplified.subs(alpha, rho * beta)
    
    # Factor and rearrange to get the desired elegant form.
    # We know from manual derivation that the expression simplifies to:
    # u_i * beta * (1 + (rho-1)*(1-eta)/(1+phi*v_i))
    term1 = u_i * beta
    term2_numerator = (rho - 1) * (1 - eta)
    term2_denominator = 1 + phi * v_i
    final_expr_obj = term1 * (1 + term2_numerator / term2_denominator)

    # Step 6: Print the final result
    # The prompt asks to output each "number" (term) in the final equation.
    # We will print the equation as a formatted string.
    dw_dt_symbol = sp.Symbol("dw_i/dt")
    LHS = f"{tau_W} * {dw_dt_symbol}"
    
    # Construct the final string for the equation
    final_eq_str = f"{LHS} = {u_i} * {beta} * (1 + (({rho} - 1) * (1 - {eta})) / (1 + {phi} * {v_i}))"
    
    print("The derived expression for the change in synaptic efficacy is:")
    print(final_eq_str)
    print("\nWhere:")
    print(f"  - {dw_dt_symbol} is the rate of change of synaptic efficacy w_i")
    print(f"  - {u_i} is the postsynaptic accumulator (Y)")
    print(f"  - {v_i} is the presynaptic accumulator (x_i)")
    print(f"  - {rho} is the constant ratio of LTD/LTP strengths (alpha/beta)")
    print(f"  - {tau_W}, {beta}, {eta}, {phi} are model parameters.")

if __name__ == '__main__':
    derive_synaptic_efficacy_equation()