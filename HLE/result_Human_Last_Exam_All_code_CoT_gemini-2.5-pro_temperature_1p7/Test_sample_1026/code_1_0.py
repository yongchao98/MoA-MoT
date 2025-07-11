import math

def calculate_guide_displacement():
    """
    Calculates the horizontal displacement of the guide.
    """
    # Given parameters
    m = 0.20  # kg, mass of the body
    M = 0.80  # kg, mass of the guide
    R_cm = 20   # cm, radius of the circular arcs
    d_cm = 50   # cm, length of the straight section
    mu_D = 0.20 # coefficient of kinetic friction

    # Convert units from cm to m
    R = R_cm / 100.0
    d = d_cm / 100.0

    # Step 1: Calculate the final height h
    # h = R - mu_D * d
    h = R - mu_D * d
    
    # Step 2: Calculate the horizontal displacement of the mass relative to the guide
    # delta_x_rel = R + d + sqrt(h * (2*R - h))
    # A simplified form is sqrt(R^2 - (mu_D*d)^2)
    sqrt_term_val = math.sqrt(h * (2 * R - h))
    delta_x_rel = R + d + sqrt_term_val

    # Step 3: Calculate the displacement of the guide D_G
    # D_G = - (m / (m + M)) * delta_x_rel
    guide_displacement = - (m / (m + M)) * delta_x_rel
    
    # Print the final equation with all the numbers
    # The term sqrt(h*(2R-h)) can also be written as sqrt(R^2 - (mu_D*d)^2)
    print("The final equation for the guide's displacement (D_G) is:")
    print("D_G = - (m / (m + M)) * [R + d + sqrt(R^2 - (mu_D * d)^2)]")
    print("\nPlugging in the numbers:")
    
    # Using f-string to display the equation with values
    equation_str = (
        f"D_G = - ({m:.2f} / ({m:.2f} + {M:.2f})) * "
        f"[{R:.2f} + {d:.2f} + sqrt({R:.2f}^2 - ({mu_D:.2f} * {d:.2f})^2)]"
    )
    print(equation_str)
    
    # Print the calculated result
    print(f"\nCalculated horizontal displacement of the guide: {guide_displacement:.4f} m")
    print(f"This is equal to {guide_displacement*100:.2f} cm to the left.")
    
    return guide_displacement

# Run the calculation and store the result
final_displacement = calculate_guide_displacement()

# Final Answer format
print(f"\n<<<{final_displacement}>>>")