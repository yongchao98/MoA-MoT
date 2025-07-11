import math

def calculate_field_strength_factor():
    """
    Demonstrates how option C leads to a forward shift of the center of gravity.

    Option C states: "Field strength varies inversely with apparent propagation
    time in any given reference frame".

    We model this by interpreting "apparent propagation time" as being modified
    by the relativistic Doppler factor (1 - n·v/c), where 'n' is the direction
    to the observer and 'v' is the source velocity. This is the core of the
    Liénard-Wiechert potential, which describes the field of a moving source.

    A field strength inversely proportional to this "time" will be stronger in
    the forward direction, indicating an effective shift of the center of gravity.

    The code calculates a factor proportional to the field strength at two points
    equidistant from a moving mass: one point in front and one behind.
    """

    # --- Setup of the scenario ---
    # Speed of light
    c = 3.0e8  # m/s
    # Velocity of the moving mass (Mass 2)
    # Let's use 80% of the speed of light
    v = 0.8 * c
    # Distance from the mass to the test points
    distance = 1.0  # meters

    beta = v / c

    # --- Calculation for the point IN FRONT of the moving mass ---
    # Here, the direction vector 'n' is parallel to the velocity 'v'.
    # So, the dot product n·(v/c) is just v/c.
    # The modifying factor is (1 - v/c).
    denominator_front = distance * (1 - beta)
    strength_factor_front = 1 / denominator_front

    # --- Calculation for the point BEHIND the moving mass ---
    # Here, the direction vector 'n' is anti-parallel to the velocity 'v'.
    # So, the dot product n·(v/c) is -v/c.
    # The modifying factor is (1 - (-v/c)) = (1 + v/c).
    denominator_back = distance * (1 + beta)
    strength_factor_back = 1 / denominator_back

    # --- Print Results ---
    print("This script demonstrates that assuming 'Field strength varies inversely with apparent propagation time'")
    print("leads to a stronger field in the direction of motion (a forward shift).\n")
    print("Parameters:")
    print(f"  Velocity (v): {v:.2e} m/s")
    print(f"  Speed of Light (c): {c:.2e} m/s")
    print(f"  Distance (d): {distance:.1f} m")
    print(f"  Beta (v/c): {beta:.1f}\n")

    print("--- Field Strength Factor (in front of mass) ---")
    print(f"Proportional to: 1 / (d * (1 - v/c))")
    print(f"Calculation: 1 / ({distance:.1f} * (1 - {beta:.1f})) = 1 / ({distance:.1f} * {1 - beta:.1f})")
    print(f"Result: {strength_factor_front:.3f}\n")

    print("--- Field Strength Factor (behind mass) ---")
    print(f"Proportional to: 1 / (d * (1 + v/c))")
    print(f"Calculation: 1 / ({distance:.1f} * (1 + {beta:.1f})) = 1 / ({distance:.1f} * {1 + beta:.1f})")
    print(f"Result: {strength_factor_back:.3f}\n")
    
    print("Conclusion: The field is significantly stronger in front of the moving mass,")
    print("which corresponds to its center of gravity appearing shifted in the direction of motion.")

calculate_field_strength_factor()
<<<C>>>