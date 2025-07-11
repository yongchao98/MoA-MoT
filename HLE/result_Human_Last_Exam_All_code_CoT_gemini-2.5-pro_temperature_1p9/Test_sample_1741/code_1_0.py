import math

def calculate_proton_position():
    """
    Calculates the final position of a proton in crossed electric and
    magnetic fields after three complete revolutions.
    """

    # --- Given constants ---
    # Charge of a proton in Coulombs
    Q_p = 1.602e-19
    # Mass of a proton in kg
    m_p = 1.673e-27
    # Magnetic field component in the z-direction in Tesla
    B_z = 5.0e-4
    # Electric field component in the x-direction in V/m
    E_x = 5.0
    # Number of revolutions
    n = 3

    # --- Calculation ---
    # The final position is determined by the total drift over n revolutions.
    # The motion along the x and z axes has no net displacement after integer revolutions,
    # and the drift velocity has no components in these directions.
    # Therefore, the final x and z positions are 0.
    final_pos_x = 0.0
    final_pos_z = 0.0

    # The final y-position is given by the formula:
    # y = (Drift Velocity_y) * (Total Time)
    # y = (-E_x / B_z) * (n * Period)
    # y = (-E_x / B_z) * (n * 2 * pi * m_p / (Q_p * B_z))
    # Simplifying this gives:
    # y = -(n * 2 * pi * m_p * E_x) / (Q_p * B_z^2)

    # --- Output the Explanation and Calculation ---
    print("The final position of the proton after three complete revolutions is determined by the total displacement due to the E-cross-B drift.")
    print("\nThe equation for the final y-position is derived from y = (drift velocity) * (total time):")
    print("y = -(n * 2 * pi * m_p * E_x) / (Q_p * B_z^2)")
    print("\nSubstituting the given values into the equation:")
    # The following line prints the equation with all the numbers substituted in, as requested.
    print(f"y = -({n} * 2 * pi * {m_p} * {E_x}) / ({Q_p} * ({B_z})**2)")

    # Calculate the final y-position using the formula
    final_pos_y = -(n * 2 * math.pi * m_p * E_x) / (Q_p * B_z**2)

    print("\nResult:")
    print(f"The final position vector of the proton is:")
    print(f"(x, y, z) = ({final_pos_x:.2f} m, {final_pos_y:.2f} m, {final_pos_z:.2f} m)")


if __name__ == '__main__':
    calculate_proton_position()