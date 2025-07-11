def solve_mystical_challenge():
    """
    Calculates the force F required to lift a rope according to the problem description.
    """

    # --- Parameters of the magical rope ---
    # You can change these values for your specific problem.
    rope_mass_m = 10.0      # mass in kilograms (kg)
    rope_length_l = 20.0    # length in meters (m)
    final_speed_v = 4.0     # speed in meters per second (m/s)

    # Standard acceleration due to gravity
    gravity_g = 9.8         # in meters per second squared (m/s^2)

    # Check for invalid input to avoid division by zero
    if rope_length_l <= 0:
        print("Error: Rope length must be a positive number.")
        return

    # Using the derived formula from the Work-Energy Theorem:
    # F = m * (g/2 + v^2/(2*l))
    force_F = rope_mass_m * (gravity_g / 2 + final_speed_v**2 / (2 * rope_length_l))

    # Output the final equation with all numbers and the result
    print(
        f"{rope_mass_m} * ({gravity_g} / 2 + {final_speed_v}**2 / (2 * {rope_length_l})) = {force_F:.2f}"
    )

solve_mystical_challenge()