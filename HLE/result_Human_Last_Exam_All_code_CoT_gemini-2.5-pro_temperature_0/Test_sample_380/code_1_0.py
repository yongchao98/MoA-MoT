import numpy as np

def calculate_field_aligned_components():
    """
    Demonstrates the transformation of magnetic field data from a standard
    coordinate system (GSE) to a Field-Aligned Coordinate (FAC) system.
    This is the correct preparatory step for calculating magnetic helicity
    of waves like AIC waves.
    """
    # --- 1. Sample Magnetic Field Data in GSE coordinates (Bx, By, Bz) in nT ---
    # This data simulates a background Parker Spiral field (~45 degrees)
    # with a superimposed wave fluctuation.
    # B_background = [-3.5, 3.5, 0.5] nT
    # Wave is a simple sine wave on the Y and Z components.
    time = np.linspace(0, 10, 100)
    # Bx(t), By(t), Bz(t)
    b_gse = np.array([
        -3.5 + 0.1 * np.sin(time),
        3.5 + 0.5 * np.sin(time * 2),
        0.5 + 0.5 * np.cos(time * 2)
    ]).T # Transpose to get a list of vectors [Bx, By, Bz] for each time step

    print("--- Step 1: Original Data ---")
    print(f"Sample magnetic field vector at t=0 in GSE (nT): {b_gse[0]}")
    print("GSE coordinates assume X is radial (Sun-Earth line), Y and Z are perpendicular.\n")

    # --- 2. Calculate the average background magnetic field, B0 ---
    b0_gse = np.mean(b_gse, axis=0)
    
    print("--- Step 2: Determine Local Background Field (B0) ---")
    print(f"Average background field B0 in GSE (nT): {np.round(b0_gse, 2)}")
    b0_magnitude = np.linalg.norm(b0_gse)
    # Angle with the radial (-X) direction
    radial_dir = np.array([-1, 0, 0])
    angle = np.arccos(np.dot(b0_gse, radial_dir) / b0_magnitude) * 180 / np.pi
    print(f"This B0 vector is at an angle of {np.round(angle, 1)} degrees to the radial direction.")
    print("This demonstrates the field is NOT radial.\n")

    # --- 3. Define the Field-Aligned Coordinate (FAC) system ---
    # z_fac is parallel to the background magnetic field B0
    z_fac = b0_gse / b0_magnitude
    # y_fac is perpendicular to the plane containing B0 and the GSE Z-axis
    # (This is a common choice to define the new system)
    y_fac = np.cross(z_fac, np.array([0, 0, 1]))
    y_fac = y_fac / np.linalg.norm(y_fac)
    # x_fac completes the right-handed system (x = y cross z)
    x_fac = np.cross(y_fac, z_fac)

    print("--- Step 3: Create Field-Aligned Coordinate (FAC) System ---")
    print(f"Parallel direction (z_fac):     {np.round(z_fac, 3)}")
    print(f"Perpendicular direction 1 (x_fac): {np.round(x_fac, 3)}")
    print(f"Perpendicular direction 2 (y_fac): {np.round(y_fac, 3)}\n")

    # --- 4. Transform the fluctuating field into the FAC system ---
    # First, get the fluctuating part of the field
    db_gse = b_gse - b0_gse
    
    # Project the fluctuations onto the new basis vectors
    # np.dot works on the whole array of vectors
    db_parallel = np.dot(db_gse, z_fac)
    db_perp1 = np.dot(db_gse, x_fac)
    db_perp2 = np.dot(db_gse, y_fac)

    print("--- Step 4: Transform Fluctuations into FAC System ---")
    print("The magnetic helicity of AIC waves should be calculated using the two PERPENDICULAR components.")
    print("Let's check the components for the first time step (t=0):")
    print(f"Original GSE fluctuations (d_Bx, d_By, d_Bz): {np.round(db_gse[0], 3)}")
    print(f"Transformed FAC components:")
    print(f"  - Parallel fluctuation (should be small for AIC): {np.round(db_parallel[0], 3)}")
    print(f"  - Perpendicular fluctuation 1 (CORRECT TO USE):   {np.round(db_perp1[0], 3)}")
    print(f"  - Perpendicular fluctuation 2 (CORRECT TO USE):   {np.round(db_perp2[0], 3)}")
    print("\nCONCLUSION: One must use the transformed perpendicular components (db_perp1, db_perp2) for helicity analysis, not the original GSE Y and Z components.")

# Run the demonstration
calculate_field_aligned_components()