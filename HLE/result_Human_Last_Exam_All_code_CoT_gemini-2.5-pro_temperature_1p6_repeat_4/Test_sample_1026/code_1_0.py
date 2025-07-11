import math

def calculate_guide_displacement():
    """
    Calculates the horizontal displacement of the guide based on the conservation of the center of mass.
    """
    # Given parameters
    m = 0.20  # kg, mass of the body
    M = 0.80  # kg, mass of the guide
    R_cm = 20   # cm, radius of the circular arcs
    d_cm = 50   # cm, length of the straight section

    # Convert units from cm to meters for calculation
    R = R_cm / 100.0
    d = d_cm / 100.0

    # The horizontal displacement of mass 'm' relative to the guide
    # is the horizontal distance from the start point to the end of the straight section.
    delta_x_mass_relative = R + d

    # Using the conservation of the center of mass principle:
    # (m + M) * delta_X_guide + m * delta_x_mass_relative = 0
    # We solve for delta_X_guide.
    delta_X_guide = - (m * delta_x_mass_relative) / (m + M)

    # Print the final equation with all the numbers
    print("The formula for the guide's displacement (ΔX_M) is derived from the conservation of the Center of Mass:")
    print("ΔX_M = - (m * (R + d)) / (m + M)")
    print("\nSubstituting the values:")
    print(f"ΔX_M = - ({m} kg * ({R} m + {d} m)) / ({m} kg + {M} kg)")
    
    # Perform the intermediate calculation steps for clarity
    numerator = - (m * (R + d))
    denominator = m + M
    
    print(f"ΔX_M = {numerator:.2f} kg*m / {denominator:.2f} kg")

    # Print the final result
    print(f"\nThe calculated horizontal displacement of the guide is {delta_X_guide:.2f} meters.")
    print(f"In centimeters, this is {delta_X_guide * 100:.0f} cm.")

calculate_guide_displacement()