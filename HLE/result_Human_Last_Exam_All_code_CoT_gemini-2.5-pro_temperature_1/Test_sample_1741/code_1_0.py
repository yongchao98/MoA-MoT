import math

def calculate_proton_position():
    """
    Calculates the final position of a proton in crossed electric and
    magnetic fields after three complete revolutions.
    """
    # --- Given values ---
    B_z = 5.0e-4      # Magnetic field in z-direction (T)
    E_x = 5.0         # Electric field in x-direction (V/m)
    Q_p = 1.602e-19   # Proton charge (C)
    m_p = 1.673e-27   # Proton mass (kg)
    N = 3             # Number of revolutions

    # --- Calculation ---
    # The final position is the result of the drift velocity multiplied by the total time.
    # The displacement from the circular part of the motion is zero after an integer number of revolutions.
    #
    # Final Position (r) = Drift Velocity (U_d) * Total Time (t_final)
    # r = ( (E x B) / B^2 ) * (N * Period)
    # r = ( (E x B) / B^2 ) * (N * 2 * pi * m / (q * B))
    #
    # For E in x-dir and B in z-dir, E x B is in -y dir.
    # So the only non-zero component is y.
    # y_final = (-E_x / B_z) * (N * 2 * pi * m_p / (Q_p * B_z))
    # y_final = - (2 * pi * N * m_p * E_x) / (Q_p * B_z^2)

    y_final_numerator = -2 * math.pi * N * m_p * E_x
    y_final_denominator = Q_p * (B_z**2)
    y_final = y_final_numerator / y_final_denominator

    # The x and z components of the final position are zero because the drift
    # is purely in the y-direction.
    x_final = 0.0
    z_final = 0.0

    # --- Output the results ---
    print("The final position of the proton can be determined by the formula:")
    print("r_final = U_drift * t_total")
    print("This simplifies to a displacement only in the y-direction:")
    print("y_final = - (2 * pi * N * m_p * E_x) / (Q_p * B_z^2)\n")

    print("Substituting the given values into the equation:")
    print(f"y_final = - (2 * {math.pi:.5f} * {N} * {m_p:.4e} kg * {E_x:.1f} V/m) / ({Q_p:.4e} C * ({B_z:.1e} T)^2)")
    print(f"y_final = ({y_final_numerator:.4e}) / ({y_final_denominator:.4e})")
    print(f"y_final = {y_final:.4f} m\n")

    print("Final Position Vector:")
    print(f"r = ({x_final:.1f} ax + {y_final:.2f} ay + {z_final:.1f} az) m")

# Execute the function
calculate_proton_position()