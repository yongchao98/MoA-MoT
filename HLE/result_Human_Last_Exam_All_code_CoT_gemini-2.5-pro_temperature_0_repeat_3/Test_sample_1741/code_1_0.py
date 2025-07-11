import math

def solve_proton_motion():
    """
    Calculates the final position of a proton in crossed electric and magnetic fields
    after three complete revolutions.
    """
    # --- Define the physical constants and initial conditions ---
    B_z = 5.0e-4      # Magnetic field in T (in the z-direction)
    E_x = 5.0         # Electric field in V/m (in the x-direction)
    Q_p = 1.602e-19   # Proton charge in C
    m_p = 1.673e-27   # Proton mass in kg
    N_rev = 3         # Number of revolutions

    # --- Explanation ---
    print("The motion of the proton is a superposition of a circular motion and a linear drift.")
    print("After an integer number of revolutions, the net displacement from the circular motion is zero.")
    print("Thus, the final position is determined by the drift velocity over the total time.\n")

    # --- Step 1: Calculate the period of one revolution (T) ---
    # The period of cyclotron motion is given by T = 2 * pi * m / (|q| * B)
    T = (2 * math.pi * m_p) / (Q_p * B_z)
    print("Step 1: Calculate the period of one revolution (T)")
    print(f"T = (2 * pi * m_p) / (Q_p * B_z)")
    print(f"T = (2 * {math.pi:.5f} * {m_p:.4e}) / ({Q_p:.4e} * {B_z:.1e}) = {T:.4e} s\n")

    # --- Step 2: Calculate the total time for three revolutions ---
    t_total = N_rev * T
    print("Step 2: Calculate the total time for three revolutions (t_total)")
    print(f"t_total = {N_rev} * T = {N_rev} * {T:.4e} s = {t_total:.4e} s\n")

    # --- Step 3: Calculate the drift velocity (U_d) ---
    # The drift velocity is U_d = (E x B) / B^2
    # U_d = (E_x * a_x) x (B_z * a_z) / B_z^2 = -(E_x / B_z) * a_y
    U_dx = 0.0
    U_dy = -E_x / B_z
    U_dz = 0.0
    print("Step 3: Calculate the drift velocity (U_d)")
    print(f"U_d = -(E_x / B_z) * a_y")
    print(f"U_d = -({E_x} / {B_z:.1e}) * a_y = {U_dy:.1e} a_y m/s\n")

    # --- Step 4: Calculate the final position (r_final) ---
    # The final position is r_final = U_d * t_total, since the proton starts at the origin.
    # The final position vector is (x_final, y_final, z_final).
    x_final = U_dx * t_total
    y_final = U_dy * t_total
    z_final = U_dz * t_total
    
    print("Step 4: Calculate the final position (r_final = U_d * t_total)")
    print("The final position has components (x, y, z):")
    print(f"x_final = {U_dx:.1f} m/s * {t_total:.4e} s = {x_final:.4f} m")
    
    # Show the final equation for the y-component with all numbers
    print(f"y_final = U_dy * t_total = (-E_x / B_z) * (N_rev * 2 * pi * m_p / (Q_p * B_z))")
    print(f"y_final = -({N_rev} * 2 * pi * m_p * E_x) / (Q_p * B_z^2)")
    print(f"y_final = -({N_rev} * 2 * {math.pi:.5f} * {m_p:.4e} * {E_x}) / ({Q_p:.4e} * ({B_z:.1e})^2)")
    print(f"y_final = {y_final:.4f} m")
    
    print(f"z_final = {U_dz:.1f} m/s * {t_total:.4e} s = {z_final:.4f} m\n")

    # --- Final Answer ---
    print("--- Final Answer ---")
    print(f"The final position of the proton is ({x_final:.4f}, {y_final:.4f}, {z_final:.4f}) m.")

solve_proton_motion()