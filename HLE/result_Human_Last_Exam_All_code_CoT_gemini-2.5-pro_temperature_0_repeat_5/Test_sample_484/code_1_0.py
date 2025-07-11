import sympy

def derive_simplified_model():
    """
    Performs a steady-state analysis of the biophysical model to derive
    a simplified expression for synaptic weight dynamics.
    """
    # 1. Define all symbolic variables from the model description.
    # We use sympy.Symbol for all variables and parameters.
    # For the derivative, we'll just use a symbol representing dot(w_i).
    tau_w, w_i = sympy.symbols('tau_w w_i')
    dot_w_i = sympy.Symbol(r'\dot{w_i}') # Represents dw_i/dt
    alpha, beta, eta = sympy.symbols('alpha beta eta')
    P_i, B_i = sympy.symbols('P_i B_i')
    M_i, Y = sympy.symbols('M_i Y')

    print("Step 1: Start with the equation for synaptic efficacy dynamics.")
    # The equation is tau_W * dW_i/dt = alpha * P_i + beta * B_i
    # We use w_i instead of W_i for simplicity.
    lhs = tau_w * dot_w_i
    rhs_initial = alpha * P_i + beta * B_i
    print(f"Original equation: {sympy.Eq(lhs, rhs_initial)}\n")

    print("Step 2: Use the steady-state expressions for P_i and B_i.")
    # From the steady-state analysis (d/dt = 0) of the equations for P_i and B_i, we get:
    # P_i = (1 - eta) * Y / (1 + M_i)
    # B_i = (eta + M_i) * Y / (1 + M_i)
    P_i_expr = (1 - eta) * Y / (1 + M_i)
    B_i_expr = (eta + M_i) * Y / (1 + M_i)
    print(f"Steady-state P_i = {P_i_expr}")
    print(f"Steady-state B_i = {B_i_expr}\n")

    print("Step 3: Substitute these expressions into the efficacy equation.")
    rhs_substituted = rhs_initial.subs({P_i: P_i_expr, B_i: B_i_expr})
    print(f"Substituted equation: {sympy.Eq(lhs, rhs_substituted)}\n")

    print("Step 4: Algebraically simplify the right-hand side.")
    # We can factor out Y / (1 + M_i)
    rhs_simplified = sympy.simplify(rhs_substituted)
    # The simplified expression is Y*(alpha*(1 - eta) + beta*(eta + M_i))/(1 + M_i)
    print(f"Simplified equation: {sympy.Eq(lhs, rhs_simplified)}\n")

    print("Step 5: Introduce the new variables u_i, v_i, and rho as defined in the problem.")
    # The problem defines:
    # v_i: presynaptic accumulator. From steady-state of M_i, M_i = phi * x_i. We set v_i = M_i.
    # u_i: postsynaptic accumulator. From steady-state of Y, Y = sum(w_j * x_j). We set u_i = Y.
    # rho: a constant. We define it to group terms.
    u_i, v_i, rho = sympy.symbols('u_i v_i rho')
    
    # Definition of rho from the expression: rho = alpha*(1-eta) + beta*eta
    rho_expr = alpha * (1 - eta) + beta * eta
    print(f"Defining presynaptic accumulator: v_i = M_i")
    print(f"Defining postsynaptic accumulator: u_i = Y")
    print(f"Defining the constant rho = {rho_expr}\n")

    print("Step 6: Rewrite the equation in terms of u_i, v_i, and rho.")
    # We substitute M_i -> v_i, Y -> u_i, and the expression for rho.
    # rhs_simplified = Y * (alpha - alpha*eta + beta*eta + beta*M_i) / (1 + M_i)
    # rhs_simplified = Y * ( (alpha*(1-eta) + beta*eta) + beta*M_i ) / (1 + M_i)
    # This becomes:
    final_rhs = u_i * (rho + beta * v_i) / (1 + v_i)
    final_equation = sympy.Eq(lhs, final_rhs)

    print("Final derived expression for the synaptic efficacy dynamics:")
    # The instruction "output each number in the final equation" is interpreted as
    # printing the final symbolic equation, as no numerical values were provided.
    print(final_equation)

if __name__ == '__main__':
    derive_simplified_model()