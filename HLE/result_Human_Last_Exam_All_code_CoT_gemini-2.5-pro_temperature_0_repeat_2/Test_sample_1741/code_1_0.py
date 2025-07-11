import math

def solve_proton_position():
    """
    Calculates the position of a proton after three revolutions in crossed
    electric and magnetic fields.
    """
    # --- Given constants ---
    Q_p = 1.602e-19  # Charge of a proton in Coulombs
    m_p = 1.673e-27  # Mass of a proton in kg
    B_z = 5.0e-4     # Magnetic field in Tesla (in z-direction)
    E_x = 5.0        # Electric field in V/m (in x-direction)
    revolutions = 3

    # The motion of the proton is a superposition of a circular motion and a linear drift.
    # After an integer number of revolutions, the net displacement is due to the E x B drift.
    # The drift velocity V_d = (E x B) / B^2.
    # Since E is in the x-direction and B is in the z-direction, the drift is in the -y direction.
    # V_dy = -E_x / B_z
    # The period of one revolution is T = 2 * pi * m_p / (Q_p * B_z).
    # The total time is t = 3 * T.
    # The final position is r = V_d * t.
    # The final y-position is y = V_dy * t = (-E_x / B_z) * (3 * 2 * pi * m_p / (Q_p * B_z))
    # This simplifies to: y = - (6 * pi * E_x * m_p) / (Q_p * B_z^2)

    # --- Calculation ---
    final_x = 0.0
    final_z = 0.0

    # Calculate the final y-position using the simplified formula
    numerator = -revolutions * 2 * math.pi * E_x * m_p
    denominator = Q_p * (B_z ** 2)
    final_y = numerator / denominator

    # --- Output ---
    print("The equation for the final y-position is:")
    print("y = - (n * 2 * pi * E_x * m_p) / (Q_p * B_z^2)")
    print("\nSubstituting the values for n=3 revolutions:")
    print(f"y = - ({revolutions} * 2 * {math.pi:.5f} * {E_x} * {m_p:.4e}) / ({Q_p:.4e} * ({B_z:.1e})^2)")
    
    print(f"\nCalculated numerator = {numerator:.4e}")
    print(f"Calculated denominator = {denominator:.4e}")

    print(f"\nFinal position (x, y, z) after {revolutions} revolutions:")
    print(f"({final_x:.4f}, {final_y:.4f}, {final_z:.4f}) meters")

solve_proton_position()