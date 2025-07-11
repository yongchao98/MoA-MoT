import numpy as np

def calculate_and_demonstrate_components():
    """
    Calculates the Parker Spiral angle at 1 AU and demonstrates how the
    perpendicular components of a magnetic field fluctuation differ when
    referenced against the radial direction vs. the local magnetic field.
    """
    # --- Part 1: Calculate the Parker Spiral Angle ---
    # Define typical parameters at the L1 point (approx. 1 AU)
    V_sw = 400e3  # Solar wind speed in m/s (400 km/s)
    r_AU = 1.496e11 # 1 Astronomical Unit in meters
    Omega_sun = 2.9e-6 # Solar angular rotation speed in rad/s

    # Calculate the Parker Spiral angle (psi)
    # tan(psi) = (Omega * r) / V_sw
    psi_rad = np.arctan((Omega_sun * r_AU) / V_sw)
    psi_deg = np.rad2deg(psi_rad)

    print("--- Parker Spiral Calculation ---")
    print(f"For a solar wind speed of {V_sw / 1e3:.0f} km/s at 1 AU:")
    print(f"The Parker Spiral angle is approximately {psi_deg:.2f} degrees from radial.")
    print("This shows the local magnetic field is NOT purely radial (which would be 0 degrees).\n")

    # --- Part 2: Demonstrate the impact on 'Perpendicular Components' ---
    # We will work in the RTN (Radial, Tangential, Normal) coordinate system.
    # Radial (R) is the x-axis, Tangential (T) is the y-axis, Normal (N) is the z-axis.

    # Define the unit vector for the RADIAL direction
    R_hat = np.array([1.0, 0.0, 0.0])

    # Define the unit vector for the LOCAL MAGNETIC FIELD (B0) based on the Parker angle.
    # The field has a radial component (cos) and a tangential component (sin).
    # The tangential component is negative in the standard RTN frame.
    B0_hat = np.array([np.cos(psi_rad), -np.sin(psi_rad), 0.0])

    # Let's imagine a spacecraft measures a magnetic field fluctuation (wave).
    # This is just a sample vector for demonstration.
    db_vec = np.array([0.1, 0.5, 0.3]) # Sample fluctuation in nT [R, T, N]

    print("--- Component Projection Demonstration ---")
    print(f"Original sample fluctuation vector db = [{db_vec[0]}, {db_vec[1]}, {db_vec[2]}]")
    print("-" * 35)

    # --- Calculation A: Components perpendicular to the RADIAL direction ---
    # This is the simplified method the user questioned.
    print("A) Using components perpendicular to the RADIAL direction:")
    
    # Project db onto the radial direction
    db_parallel_R_component = np.dot(db_vec, R_hat)
    # The full parallel vector is: db_parallel_R_component * R_hat
    
    # The perpendicular component is the original vector minus the parallel part.
    db_perp_R = db_vec - (db_parallel_R_component * R_hat)
    
    # In this frame, the "perpendicular components" are the T and N components.
    print(f"Component parallel to Radial: {db_parallel_R_component:.4f}")
    print(f"Resulting perpendicular vector db_perp_R = [{db_perp_R[0]:.4f}, {db_perp_R[1]:.4f}, {db_perp_R[2]:.4f}]")
    print("The helicity calculation would use the components [0.5000, 0.3000].\n")


    # --- Calculation B: Components perpendicular to the LOCAL MAGNETIC FIELD ---
    # This is the more physically rigorous method.
    print("B) Using components perpendicular to the LOCAL MAGNETIC FIELD (B0):")
    
    # Project db onto the B0 direction
    db_parallel_B0_component = np.dot(db_vec, B0_hat)
    
    # The perpendicular component is the original vector minus the parallel part.
    db_perp_B0 = db_vec - (db_parallel_B0_component * B0_hat)
    
    print(f"Component parallel to B0: {db_parallel_B0_component:.4f}")
    print(f"Resulting perpendicular vector db_perp_B0 = [{db_perp_B0[0]:.4f}, {db_perp_B0[1]:.4f}, {db_perp_B0[2]:.4f}]")
    print("The helicity calculation would use the components of this new vector.")
    print("-" * 35)

    print("\nCONCLUSION:")
    print("As shown, the vector representing the 'perpendicular components' is different in each case.")
    print("This confirms your intuition that using a radial frame is an approximation, and a physically rigorous analysis should be done relative to the local magnetic field.")

if __name__ == '__main__':
    calculate_and_demonstrate_components()