import math

def calculate_proton_position():
    """
    Calculates the position of a proton after three complete revolutions in
    crossed electric and magnetic fields.
    """
    # --- Given constants ---
    B_z = 5.0e-4      # Magnetic field in z-direction (T)
    E_x = 5.0         # Electric field in x-direction (V/m)
    Q_p = 1.602e-19   # Proton charge (C)
    m_p = 1.673e-27   # Proton mass (kg)
    # The initial velocity does not affect the final position after an integer
    # number of revolutions, so it's not needed for the final calculation.
    # U0_x = 2.5e5 m/s

    # --- Calculations ---
    # The motion is a cycloid. After an integer number of revolutions, the
    # displacement is due to the E x B drift.
    # The final position (x, y, z) is determined as follows:
    # x(final) = 0
    # z(final) = 0
    # y(final) = drift_velocity_y * total_time
    #
    # Drift velocity (v_dy) = -E_x / B_z
    # Angular frequency (omega_c) = (Q_p * B_z) / m_p
    # Period of one revolution (T) = 2 * pi / omega_c
    # Total time for 3 revolutions (t) = 3 * T = 6 * pi / omega_c
    #
    # y(final) = (-E_x / B_z) * (6 * pi / omega_c)
    # Substituting omega_c:
    # y(final) = (-E_x / B_z) * (6 * pi * m_p) / (Q_p * B_z)
    # y(final) = -(6 * pi * m_p * E_x) / (Q_p * B_z^2)

    # Final position components
    x_final = 0.0
    z_final = 0.0
    
    # Calculate the final y-position
    numerator = 6 * math.pi * m_p * E_x
    denominator = Q_p * (B_z ** 2)
    y_final = -numerator / denominator

    # --- Output the results ---
    print("The final position of the proton is determined by the drift over three revolution periods.")
    print("The formula for the final y-position is: y = -(6 * pi * m_p * E_x) / (Q_p * B_z^2)\n")
    
    print("Plugging in the numbers:")
    print(f"y = - (6 * {math.pi:.4f} * {m_p:.4e} * {E_x:.1f}) / ({Q_p:.4e} * ({B_z:.1e})^2)")
    
    # Calculate and show the intermediate values for clarity
    print(f"y = - ({numerator:.4e}) / ({denominator:.4e})")

    print(f"\nFinal position after three complete revolutions:")
    print(f"x = {x_final:.2f} m")
    print(f"y = {y_final:.2f} m")
    print(f"z = {z_final:.2f} m")
    
    # Final answer in the required format
    final_pos_str = f"({x_final:.2f}, {y_final:.2f}, {z_final:.2f}) m"
    print(f"\n<<< {final_pos_str} >>>")


calculate_proton_position()