import math

def calculate_guide_displacement():
    """
    Calculates the horizontal displacement of the guide.
    """
    # Given parameters
    m = 0.20  # kg, mass of the body
    M = 0.80  # kg, mass of the guide
    R = 0.20  # m, radius of the circular arcs (20 cm)
    d = 0.50  # m, length of the straight section (50 cm)
    mu_D = 0.20 # coefficient of dynamic friction

    # Step 1: Calculate the final height 'h' using the work-energy theorem.
    # The loss in potential energy equals the work done by friction.
    # mg(R - h) = mu_D * mg * d
    # R - h = mu_D * d
    # h = R - mu_D * d
    h = R - mu_D * d

    # Step 2: Calculate the horizontal displacement of the mass 'm' relative to the guide.
    # This is the sum of the horizontal projections of the path segments:
    # 1. First arc: R
    # 2. Straight section: d
    # 3. Second arc (up to height h): sqrt(2*R*h - h^2)
    # Total relative displacement: delta_x_m_rel = R + d + sqrt(2*R*h - h^2)
    horizontal_travel_on_second_arc = math.sqrt(2 * R * h - h**2)
    delta_x_m_rel = R + d + horizontal_travel_on_second_arc

    # Step 3: Calculate the horizontal displacement of the guide 'M'.
    # From conservation of the center of mass horizontal position:
    # delta_X_M = - (m / (m + M)) * delta_x_m_rel
    delta_X_M = - (m / (m + M)) * delta_x_m_rel

    # Print the equation with all the numerical values
    print("The final height 'h' reached by the mass is calculated as:")
    print(f"h = R - mu_D * d = {R:.2f} - {mu_D:.2f} * {d:.2f} = {h:.2f} m\n")

    print("The displacement of the guide is calculated using the formula:")
    print("Displacement = - (m / (m + M)) * (R + d + sqrt(2*R*h - h^2))")
    print("\nSubstituting the numerical values into the equation:")
    print(f"Displacement = - ({m:.2f} / ({m:.2f} + {M:.2f})) * ({R:.2f} + {d:.2f} + sqrt(2*{R:.2f}*{h:.2f} - {h:.2f}**2))")

    # Print the final calculated result
    print(f"\nFinal Result:")
    print(f"The horizontal displacement of the guide is {delta_X_M:.4f} meters.")

calculate_guide_displacement()