import math

def calculate_proton_position():
    """
    Calculates the final position of a proton in crossed E and B fields
    after three complete revolutions.
    """
    # Given values and physical constants
    Qp = 1.602e-19  # Proton charge in Coulombs
    mp = 1.673e-27  # Proton mass in kg
    Bz = 5.0e-4     # Magnetic field in Tesla (in the +z direction)
    Ex = 5.0        # Electric field in V/m (in the +x direction)
    N_revolutions = 3

    print("This script calculates the final position of the proton.")
    print("-" * 50)

    # Step 1: Calculate the drift velocity, Ud = (E x B) / B^2.
    # E is along ax, B is along az, so E x B is along -ay.
    # The only non-zero component is in the y-direction.
    U_drift_y = -Ex / Bz
    print("Step 1: Calculate the drift velocity (Ud).")
    print(f"The equation for the drift velocity component is Ud_y = -E_x / B_z")
    print(f"Ud_y = -({Ex}) / ({Bz}) = {U_drift_y} m/s")
    print("-" * 50)

    # Step 2: Calculate the period of one revolution (T).
    # The angular frequency (cyclotron frequency) is omega_c = |q|*B/m.
    omega_c = Qp * Bz / mp
    # The period is T = 2*pi / omega_c.
    T_period = (2 * math.pi) / omega_c
    print("Step 2: Calculate the period of one revolution (T).")
    print(f"The equation for the period is T = 2 * pi * m / (q * B)")
    print(f"T = (2 * {math.pi:.4f} * {mp:.3e}) / ({Qp:.3e} * {Bz}) = {T_period:.4e} s")
    print("-" * 50)
    
    # Step 3: Calculate the total time for N revolutions.
    t_final = N_revolutions * T_period
    print("Step 3: Calculate the total time for 3 revolutions (t_final).")
    print(f"The equation is t_final = N * T")
    print(f"t_final = {N_revolutions} * {T_period:.4e} s = {t_final:.4e} s")
    print("-" * 50)

    # Step 4: Calculate the final position.
    # After an integer number of revolutions, the gyration displacement is zero.
    # The net displacement is from the drift: r_final = Ud * t_final.
    x_final = 0.0
    y_final = U_drift_y * t_final
    z_final = 0.0
    print("Step 4: Calculate the final position (r_final).")
    print("After 3 full revolutions, the displacement due to the circular motion is zero.")
    print("The final position is determined by the drift velocity and total time.")
    print("r_final_x = 0.0 m")
    print(f"r_final_y = Ud_y * t_final = {U_drift_y} m/s * {t_final:.4e} s = {y_final:.4f} m")
    print("r_final_z = 0.0 m")
    print("-" * 50)

    print("Final Result:")
    print(f"The final position of the proton is ({x_final:.1f}, {y_final:.4f}, {z_final:.1f}) m.")


calculate_proton_position()