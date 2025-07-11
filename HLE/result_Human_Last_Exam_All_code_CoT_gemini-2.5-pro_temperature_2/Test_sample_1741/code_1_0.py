import math

def calculate_proton_position():
    """
    Calculates the final position of a proton in crossed E and B fields
    after three complete revolutions.
    """

    # --- Given Constants ---
    Q_p = 1.602e-19  # Proton charge in Coulombs
    m_p = 1.673e-27  # Proton mass in kg
    B_z = 5.0e-4     # Magnetic field in the z-direction in Tesla
    E_x = 5.0        # Electric field in the x-direction in V/m
    revolutions = 3

    # The motion consists of a circular part and a drift part. After complete
    # revolutions, the displacement from the circular motion is zero. The final
    # position is determined by the drift velocity and the total time.

    # --- Step 1: Analyze motion in x and z directions ---
    # The E x B drift is purely in the y-direction. The circular motion is centered
    # on the drift path. After full revolutions, the particle returns to the drift path.
    # No forces act in the z-direction.
    x_final = 0.0
    z_final = 0.0

    print("The final position of the proton is given by the vector r = (x, y, z).")

    print("\n--- Final x-position ---")
    print("After three complete revolutions, the displacement from the circular component of motion in the x-direction is zero.")
    print(f"x = {x_final:.1f} m")

    print("\n--- Final z-position ---")
    print("There are no forces or initial velocity in the z-direction.")
    print(f"z = {z_final:.1f} m")

    # --- Step 2: Calculate the final position in the y-direction ---
    # The final y-position is due to the E x B drift over the total time.
    # y_final = drift_velocity_y * total_time
    # drift_velocity_y = -E_x / B_z
    # total_time = revolutions * Period = revolutions * (2 * pi / omega)
    # omega = Q_p * B_z / m_p
    # So, y_final = (-E_x / B_z) * (revolutions * 2 * pi * m_p / (Q_p * B_z))

    y_final = -(revolutions * 2 * math.pi * E_x * m_p) / (Q_p * B_z**2)
    numerator = -(revolutions * 2 * math.pi * E_x * m_p)
    denominator = Q_p * B_z**2

    print("\n--- Final y-position ---")
    print("The final y-position is calculated using the formula:")
    print("y = -(N_revs * 2 * pi * E_x * m_p) / (Q_p * B_z^2)")
    print("\nSubstituting the values:")
    print(f"y = -({revolutions} * 2 * {math.pi:.5f} * {E_x:.1f} * {m_p:.4e}) / ({Q_p:.4e} * ({B_z:.1e})^2)")
    print(f"y = ({numerator:.5e}) / ({denominator:.5e})")
    print(f"y = {y_final:.4f} m")

    print("\n-------------------------------------------")
    print("            FINAL RESULT")
    print("-------------------------------------------")
    print(f"The final position vector after {revolutions} revolutions is:")
    print(f"r = ({x_final:.1f}, {y_final:.4f}, {z_final:.1f}) m")

# Execute the function to print the solution
calculate_proton_position()