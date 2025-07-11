import numpy as np

def calculate_field_aligned_components():
    """
    This function demonstrates how to find the correct magnetic field components
    for calculating helicity by transforming into a Field-Aligned Coordinate (FAC) system.

    We start with hypothetical measurements in a spacecraft's Radial-Tangential-Normal
    (RTN) coordinate system, where 'x' is Radial.
    """
    
    # --- Step 1: Define Measured Magnetic Field Vectors ---
    # Assume we measure the average background magnetic field (B0) at L1.
    # Note it is NOT purely radial (i.e., not [B, 0, 0]). This is a typical Parker Spiral orientation.
    # Units are in nanoTeslas (nT).
    # B0 = [Br, Bt, Bn]
    B0 = np.array([-4.5, 2.5, 0.8])

    # Assume at some instant, we measure the magnetic field perturbation of the wave (delta_B).
    # delta_B = [dBr, dBt, dBn]
    delta_B = np.array([0.2, 0.7, -0.5])

    print("--- Initial Measurements (RTN coordinates) ---")
    print(f"Background Magnetic Field B0 (R,T,N): {B0} nT")
    print(f"Wave Perturbation delta_B (R,T,N):  {delta_B} nT\n")
    
    # --- Step 2: Define the Field-Aligned Coordinate (FAC) System basis vectors ---
    # The 'z_fac' axis is parallel to the local background magnetic field B0.
    z_fac = B0 / np.linalg.norm(B0)
    
    # To define the perpendicular plane, we need another vector that is not parallel to B0.
    # The radial direction is a convenient choice.
    radial_vector = np.array([1.0, 0.0, 0.0])
    
    # The 'y_fac' axis is perpendicular to the plane formed by B0 and the radial vector.
    y_fac = np.cross(z_fac, radial_vector)
    y_fac = y_fac / np.linalg.norm(y_fac)
    
    # The 'x_fac' axis completes the right-handed system (x = y x z).
    x_fac = np.cross(y_fac, z_fac)

    print("--- Field-Aligned Coordinate (FAC) System Basis Vectors ---")
    print(f"x_fac (Perpendicular 1): {np.round(x_fac, 3)}")
    print(f"y_fac (Perpendicular 2): {np.round(y_fac, 3)}")
    print(f"z_fac (Parallel to B0):    {np.round(z_fac, 3)}\n")

    # --- Step 3: Project the Wave Perturbation onto the FAC System ---
    # We use the dot product to find the components of delta_B along the new FAC axes.
    db_perp1 = np.dot(delta_B, x_fac)
    db_perp2 = np.dot(delta_B, y_fac)
    db_parallel = np.dot(delta_B, z_fac)

    print("--- Wave Components in the Correct Physical System (FAC) ---")
    print(f"Perpendicular Component 1 (along x_fac): {db_perp1:.4f} nT")
    print(f"Perpendicular Component 2 (along y_fac): {db_perp2:.4f} nT")
    print(f"Parallel Component (along z_fac):        {db_parallel:.4f} nT\n")

    # --- Step 4: Use the Correct Components for Physical Calculations ---
    # For helicity and wave power calculations, ONLY the perpendicular components are used.
    # The total power in the perpendicular plane is the sum of their squares.
    perp_power = db_perp1**2 + db_perp2**2
    
    print("--- Final Calculation Example (Perpendicular Power) ---")
    print("The helicity calculation MUST use the components perpendicular to the local magnetic field.")
    print(f"The perpendicular power is calculated from these components:")
    print(f"Perp_Power = (db_perp1)^2 + (db_perp2)^2")
    print(f"Perp_Power = ({db_perp1:.4f})^2 + ({db_perp2:.4f})^2 = {perp_power:.4f} nT^2")


if __name__ == '__main__':
    calculate_field_aligned_components()
