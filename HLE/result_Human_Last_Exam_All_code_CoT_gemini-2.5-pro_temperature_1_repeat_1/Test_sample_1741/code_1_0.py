import math

def calculate_proton_position():
    """
    Calculates the position of a proton after three revolutions in E and B fields.
    """
    # Given constants
    B_z = 5.0e-4  # Magnetic field in Tesla
    E_x = 5.0     # Electric field in V/m
    Q_p = 1.602e-19 # Proton charge in Coulombs
    m_p = 1.673e-27 # Proton mass in kg
    U_0x = 2.5e5   # Initial velocity in m/s
    revolutions = 3

    # The motion along the z-axis is zero as there are no forces in that direction.
    z_final = 0.0

    # After an integer number of revolutions, the cyclical part of the motion
    # in the xy-plane results in zero net displacement from the guiding center's path.
    # Therefore, the final x-position is also zero.
    x_final = 0.0
    
    # The net displacement is due to the E x B drift.
    # The drift velocity is v_d = - (E_x / B_z) in the y-direction.
    # The angular frequency is omega = (Q_p * B_z) / m_p.
    # The period of one revolution is T = 2 * pi / omega.
    # The time for 3 revolutions is t = 3 * T = 6 * pi / omega = (6 * pi * m_p) / (Q_p * B_z).
    # The final y position is y = v_d * t = -(E_x / B_z) * t
    # y = -(E_x / B_z) * (6 * pi * m_p) / (Q_p * B_z)
    # y = -(6 * pi * m_p * E_x) / (Q_p * B_z^2)

    y_final = -(6 * math.pi * m_p * E_x) / (Q_p * B_z**2)
    
    print("The final position of the proton is (x, y, z).")
    print(f"x = {x_final:.2f} m")
    print(f"z = {z_final:.2f} m")
    print("\nCalculation for the y-coordinate:")
    print("The formula for the final y-position is: y = -(6 * pi * m_p * E_x) / (Q_p * B_z^2)")
    print(f"y = - (6 * {math.pi:.4f} * {m_p:.4g} * {E_x:.1f}) / ({Q_p:.4g} * ({B_z:.1g})^2)")
    print(f"y = {y_final:.3f} m")
    
    print(f"\nFinal Position Vector: ({x_final:.2f}, {y_final:.2f}, {z_final:.2f}) m")

calculate_proton_position()