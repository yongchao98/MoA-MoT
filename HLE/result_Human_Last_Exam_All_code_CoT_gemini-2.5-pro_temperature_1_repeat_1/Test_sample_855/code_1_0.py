# [Precise Landing with Superconducting Computer]

def titan_calculate_force():
    """
    Simulates the force calculation on the Titan computer,
    following its 5-bit fractional arithmetic constraints.
    """
    # Step 1: Define initial values in real-world units
    probe_mass_val = 50  # kg
    initial_velocity_val = -300  # m/s
    landing_height_val = 5000  # m

    # Step 2: Calculate the required net acceleration (a_net)
    # a_net = -v0^2 / (2 * h) = -(-300)^2 / (2 * -5000) = 90000 / 10000 = 9 m/s^2
    # This calculation is simple. Representing 9 as a Titan fraction is (9, 1, 0).
    a_net_val = 9

    # Step 3: Attempt to calculate gravitational acceleration (g)
    # As explained in the plan, calculating g precisely involves intermediate
    # multiplications that result in numerators > 31 (e.g., 28 * 2 = 56).
    # This violates Titan's constraints.

    # Step 4: Apply simplification by eliminating the negligible term 'g'
    # The true value of g (~0.17 m/s^2) is negligible compared to a_net (9 m/s^2).
    # Adding g to a_net would also cause overflows (e.g., 9 + 1/6 = 55/6).
    # Therefore, the only possible calculation on Titan is to approximate g as 0.
    g_val_approx = 0

    # Step 5: Calculate the final rocket force
    # F_rocket = m * (g_approx + a_net)
    # On Titan, this is (5, 1, 1) * ((0,1,0) + (9,1,0))
    # which simplifies to (5, 1, 1) * (9, 1, 0)
    # Intermediate multiplication: (5*9, 1*1, 1+0) -> (45, 1, 1)
    # The numerator 45 exceeds 31, requiring renormalization.
    # 45 * 10^1 = 450.
    # To renormalize, we find a valid fraction and exponent: 4.5 * 10^2.
    # The fraction for 4.5 is 9/2. The exponent is 2.
    # Final Titan representation: (9, 2, 2)
    # The real value of this result is 450.
    force_val = probe_mass_val * (g_val_approx + a_net_val)

    # Step 6: Print the final equation with the calculated values
    print("Final Equation:")
    print(f"F_rocket = m_probe * (g + a_net)")
    print("\nTitan Calculation:")
    print("Due to 5-bit constraints, g must be approximated as 0.")
    print(f"F_rocket = {probe_mass_val} kg * ({g_val_approx} m/s^2 + {a_net_val} m/s^2) = {int(force_val)} N")

titan_calculate_force()