import numpy as np

def calculate_helicity_comparison():
    """
    Demonstrates the difference between calculating wave properties in a
    Field-Aligned Coordinate (FAC) system versus a fixed GSE-like system.
    """
    print("This script demonstrates why the coordinate system is crucial for helicity calculation.")
    print("-" * 70)

    # 1. Define the physical environment at L1
    # The magnetic field B is NOT radial. Due to the Parker Spiral, at 1 AU
    # it has both radial (X) and tangential (Y) components.
    # Let's assume a typical 45-degree angle in the XY-plane (ecliptic).
    # We use a GSE-like system: X=Sun->Earth, Y=duskward, Z=North.
    # Solar wind flows roughly along -X, so the B field points away from the Sun.
    # We set Bx and By components to be equal in magnitude.
    B_mean = np.array([-5.0, 5.0, 0.0])  # in nanoTesla (nT)
    B_magnitude = np.linalg.norm(B_mean)
    radial_direction = np.array([1.0, 0.0, 0.0])

    angle_deg = np.arccos(np.dot(B_mean, -radial_direction) / B_magnitude) * (180.0 / np.pi)
    print(f"Step 1: Define a realistic mean magnetic field B_mean at L1.")
    print(f"B_mean (GSE-like frame) = {B_mean} nT")
    print(f"The angle between the magnetic field and the radial direction is {angle_deg:.2f} degrees.\n")

    # 2. Create the physically correct Field-Aligned Coordinate (FAC) system
    e_parallel = B_mean / B_magnitude
    # Define a second vector to create the perpendicular plane. The ecliptic normal (GSE Z-axis) is a good choice.
    e_perp1 = np.cross(e_parallel, [0, 0, 1])
    e_perp1 /= np.linalg.norm(e_perp1)
    # The third basis vector completes the right-hand set.
    e_perp2 = np.cross(e_parallel, e_perp1)

    print("Step 2: Construct the correct Field-Aligned Coordinate (FAC) system based on B_mean.")
    print(f"  - Parallel axis (along B): {np.round(e_parallel, 3)}")
    print(f"  - Perpendicular axis 1:    {np.round(e_perp1, 3)}")
    print(f"  - Perpendicular axis 2:    {np.round(e_perp2, 3)}\n")

    # 3. Create a perfect, left-hand circularly polarized wave IN THE FAC SYSTEM
    # This is what a pure AIC wave would look like. It has no component along B_mean.
    # The components in the perpendicular plane (B_perp1, B_perp2) are equal and 90 degrees out of phase.
    db_fac = np.array([1.0, 1.0j, 0.0]) # Represents a snapshot in frequency space; e.g., cos(wt) + i*sin(wt)

    # Transform this physical wave back into the fixed (GSE) frame that a spacecraft measures.
    # This is done by projecting the FAC components onto the GSE axes.
    db_gse = db_fac[0] * e_perp1 + db_fac[1] * e_perp2 + db_fac[2] * e_parallel

    print("Step 3: Define a pure AIC wave and see its components in each frame.")

    # --- Analysis in the correct FAC frame ---
    print("\n--- In the Correct (FAC) Frame ---")
    print("The wave components are defined perpendicular to the local magnetic field B_mean.")
    print("The components used for the helicity equation are:")
    # For the equation, we show the real part of the wave snapshot
    print(f"  B_perp1 = {np.real(db_fac[0]):.4f}")
    print(f"  B_perp2 = {np.real(db_fac[1]):.4f} (Ignoring phase for this snapshot)")
    print(f"Final Equation Form: σ_m = 2 * Im(S_perp1_perp2) / (S_perp1_perp1 + S_perp2_perp2)")
    print("Result: Since the wave is perfectly circular in this frame, |σ_m| would be calculated as 1.0.")

    # --- Analysis in the simplified fixed frame ---
    print("\n--- In the Simplified (GSE) Frame ---")
    print("The calculation uses components perpendicular to the RADIAL (X) direction.")
    print("The components used for the helicity equation are:")
    print(f"  B_y = {np.real(db_gse[1]):.4f}")
    print(f"  B_z = {np.real(db_gse[2]):.4f}")
    print(f"Final Equation Form: σ_m' = 2 * Im(S_yz) / (S_yy + S_zz)")
    # Let's quantify the error.
    power_in_yz_plane = np.real(db_gse[1])**2 + np.real(db_gse[2])**2
    total_power = np.real(np.dot(db_gse, np.conj(db_gse)))
    print(f"Result: The power in the Y-Z plane is {power_in_yz_plane / total_power:.2%}, which is less than the total wave power.")
    print("This shows the fixed frame fails to capture the full wave and gives an inaccurate result.")
    print("\nConclusion: For a physically rigorous calculation of helicity, one MUST transform to a field-aligned coordinate system.")

# Run the demonstration
calculate_helicity_comparison()