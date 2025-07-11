import math

def calculate_guide_displacement():
    """
    Calculates the horizontal displacement of the guide based on the conservation of the center of mass.
    """
    # Given parameters
    m = 0.20  # kg (mass of the body)
    M = 0.80  # kg (mass of the guide)
    R = 0.20  # m (radius of the circular arcs)
    d = 0.50  # m (length of the straight section)
    mu_D = 0.20 # coefficient of dynamic friction

    # --- Step 1: Calculate the final height 'h' using the work-energy theorem ---
    # The loss in potential energy equals the work done by friction.
    # m*g*R = m*g*h + mu_D*m*g*d  =>  h = R - mu_D*d
    h = R - mu_D * d

    # --- Step 2: Calculate the horizontal displacement of the mass relative to the guide (Δx_m_rel) ---
    # We define a coordinate system with the origin at the bottom-left corner of the guide.
    # The mass starts at x_initial_rel = 0 (at the top of the first arc).
    # The mass ends on the second arc at height y = h. The center of the second arc is at (R+d, R).
    # The equation for the second arc is (x - (R+d))^2 + (y - R)^2 = R^2.
    # We solve for x when y = h:
    # x_final_rel = (R + d) - sqrt(R^2 - (h - R)^2)
    # Since h - R = -mu_D * d, this simplifies to:
    x_final_rel = (R + d) - math.sqrt(R**2 - (mu_D * d)**2)
    x_initial_rel = 0
    delta_x_m_rel = x_final_rel - x_initial_rel

    # --- Step 3: Calculate the displacement of the guide (ΔX_M) ---
    # From conservation of the center of mass: ΔX_M = - (m / (m + M)) * Δx_m_rel
    total_mass = m + M
    mass_ratio = m / total_mass
    delta_X_M = -mass_ratio * delta_x_m_rel

    # --- Output the results ---
    print("This script calculates the horizontal displacement of the guide.")
    print("\nFirst, we find the horizontal displacement of the mass relative to the guide.")
    print(f"The final relative horizontal position of the mass is {delta_x_m_rel:.4f} m.")
    
    print("\nNext, we use the conservation of the center of mass to find the guide's displacement (ΔX_M).")
    print("The formula is: ΔX_M = - (m / (m + M)) * Δx_m_rel")
    
    print("\nPlugging in the numbers:")
    print(f"ΔX_M = - ({m} / ({m} + {M})) * {delta_x_m_rel:.4f}")
    print(f"ΔX_M = - ({mass_ratio:.2f}) * {delta_x_m_rel:.4f}")
    print(f"ΔX_M = {delta_X_M:.4f} m")

    print(f"\nThe final horizontal displacement of the guide is {delta_X_M:.4f} meters, or {delta_X_M * 100:.2f} cm to the left.")
    
    # Return the final numerical answer for the specified format
    return delta_X_M

# Run the calculation and capture the final answer
final_displacement = calculate_guide_displacement()
# The final answer is requested in a specific format.
# Let's provide the value in meters, rounded to 3 significant figures.
final_answer = round(final_displacement, 3)
print(f"\n<<<__{final_answer}__>>>")
