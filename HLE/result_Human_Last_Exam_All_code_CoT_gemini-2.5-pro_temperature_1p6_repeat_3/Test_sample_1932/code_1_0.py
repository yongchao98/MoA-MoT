import math

def calculate_weight_change(d, h, H, rho, t, g=9.8):
    """
    Calculates the change in weight of an hourglass based on the dominant effect (impact force).

    The change in weight Delta_W is estimated by the impact force of the falling sand.
    Formula: Delta_W = m_dot * v_fall
    
    m_dot is the mass flow rate:
    Total sand mass M_sand = rho * Volume = rho * (pi * d^2 / 4) * h
    m_dot = M_sand / t = (pi * d^2 * h * rho) / (4 * t)
    
    v_fall is the velocity of sand upon impact.
    It falls a distance L. We approximate L as H - h.
    v_fall = sqrt(2 * g * L) = sqrt(2 * g * (H - h))
    
    Combining these gives the final expression for Delta_W.
    """
    
    # Calculate mass flow rate (m_dot)
    m_dot = (math.pi * d**2 * h * rho) / (4 * t)
    
    # Calculate fall velocity (v_fall)
    # Using L = H - h as the representative fall distance
    if H <= h:
        # This case is physically unrealistic but we handle it.
        # Fall distance would be negative, which is not possible. Let's assume zero change.
        v_fall = 0
    else:
        v_fall = math.sqrt(2 * g * (H - h))
        
    # Calculate the change in weight (Delta_W)
    delta_W = m_dot * v_fall
    
    print("This python code calculates the estimated change in weight of the hourglass.")
    print("The formula for the change in weight (Delta W) is derived from the impact force of the sand.")
    print("\nFormula components:")
    print(f"Mass flow rate (m_dot) = (pi * d^2 * h * rho) / (4 * t)")
    print(f"Fall velocity (v_fall) = sqrt(2 * g * (H - h))")
    print("\nFinal equation for Delta W:")
    print("Delta W = ( (pi * d^2 * h * rho) / (4 * t) ) * sqrt(2 * g * (H - h))")
    
    # Let's print the expression corresponding to the answer choice
    print("\nThe symbolic expression is:")
    # We will print the components of the expression from Answer A
    print("Delta W = (pi * d^2 * h * rho / (4 * t)) * sqrt(2 * g * (H - h))")

# You can run this function with the provided values, although the problem asks for the symbolic expression.
# d_val = 0.01  # m
# h_val = 0.02  # m
# H_val = 0.04  # m
# rho_val = 1500  # kg/m^3
# t_val = 60  # s
# calculate_weight_change(d_val, h_val, H_val, rho_val, t_val)

calculate_weight_change(d="d", h="h", H="H", rho="rho", t="t", g="g")