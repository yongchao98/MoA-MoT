import numpy as np

def calculate_field_aligned_components(b_data):
    """
    Transforms magnetic field data into a field-aligned coordinate system.

    The function first calculates the mean magnetic field vector to define the
    parallel direction. It then constructs a new basis (e_parallel, e_perp1, e_perp2)
    and projects the magnetic field fluctuations onto this new basis.

    Args:
        b_data (np.ndarray): An array of magnetic field vectors, shape (N, 3),
                             where N is the number of time points. Each row is [Bx, By, Bz].

    Returns:
        tuple: A tuple containing:
            - b_parallel (np.ndarray): Fluctuation component parallel to the mean field.
            - b_perp1 (np.ndarray): First perpendicular fluctuation component.
            - b_perp2 (np.ndarray): Second perpendicular fluctuation component.
    """
    # 1. Calculate the mean magnetic field B0
    b_mean = np.mean(b_data, axis=0)
    print("--- Calculating Field-Aligned Frame ---")
    print(f"Mean Magnetic Field B0 = [{b_mean[0]:.2f}, {b_mean[1]:.2f}, {b_mean[2]:.2f}] nT")

    # 2. Define the new basis vectors for the Field-Aligned Coordinate (FAC) system
    e_parallel = b_mean / np.linalg.norm(b_mean)
    
    # Create a temporary vector that is not parallel to e_parallel
    # We use the Z-axis of the original system ([0, 0, 1]) as a reference
    # If e_parallel is along Z, use Y-axis instead to avoid cross product being zero.
    ref_vec = np.array([0., 0., 1.])
    if np.abs(np.dot(e_parallel, ref_vec)) > 0.99:
        ref_vec = np.array([0., 1., 0.])

    # The first perpendicular axis is orthogonal to both B0 and the reference vector
    e_perp1 = np.cross(e_parallel, ref_vec)
    e_perp1 /= np.linalg.norm(e_perp1)

    # The second perpendicular axis completes the right-hand system
    e_perp2 = np.cross(e_parallel, e_perp1)

    print("\nNew Field-Aligned Basis Vectors:")
    print(f"  e_parallel = [{e_parallel[0]:.2f}, {e_parallel[1]:.2f}, {e_parallel[2]:.2f}]")
    print(f"  e_perp1    = [{e_perp1[0]:.2f}, {e_perp1[1]:.2f}, {e_perp1[2]:.2f}]")
    print(f"  e_perp2    = [{e_perp2[0]:.2f}, {e_perp2[1]:.2f}, {e_perp2[2]:.2f}]")

    # 3. Calculate magnetic field fluctuations
    b_fluctuations = b_data - b_mean
    print("\nExample original fluctuation vector (first data point):")
    print(f"  delta_B = [{b_fluctuations[0,0]:.2f}, {b_fluctuations[0,1]:.2f}, {b_fluctuations[0,2]:.2f}]")

    # 4. Project fluctuations onto the new basis
    # This is equivalent to a rotation: B_fac = R * B_original
    # where R is the rotation matrix with rows e_parallel, e_perp1, e_perp2
    b_parallel_comp = np.dot(b_fluctuations, e_parallel)
    b_perp1_comp = np.dot(b_fluctuations, e_perp1)
    b_perp2_comp = np.dot(b_fluctuations, e_perp2)
    
    print("\nSame fluctuation vector transformed to FAC system:")
    print(f"  delta_B_parallel = {b_parallel_comp[0]:.2f}")
    print(f"  delta_B_perp1    = {b_perp1_comp[0]:.2f}")
    print(f"  delta_B_perp2    = {b_perp2_comp[0]:.2f}")
    
    return b_parallel_comp, b_perp1_comp, b_perp2_comp


if __name__ == '__main__':
    # --- Part 1: Using the fixed coordinate system (Approximate Method) ---
    # Imagine this is magnetic field data from a spacecraft at L1 in a GSE-like frame (X is radial).
    # We create synthetic data where the mean field is tilted ~45 deg in the X-Y plane (Parker Spiral)
    time = np.linspace(0, 2 * np.pi, 100)
    # Background field components (tilted)
    b0x, b0y, b0z = 4.0, 4.0, 0.5 
    # Wave fluctuation components (circular polarization in Y-Z plane for simplicity of example)
    db_y = 1.0 * np.sin(time)
    db_z = 1.0 * np.cos(time)
    db_x = 0.1 * np.sin(time) # small parallel fluctuation
    
    b_gse = np.vstack([b0x + db_x, b0y + db_y, b0z + db_z]).T

    print("--- Part 1: Analysis in Fixed Frame (e.g., GSE/RTN) ---")
    print("This method assumes the perpendicular components are By and Bz.")
    b_gse_fluctuations = b_gse - np.mean(b_gse, axis=0)
    by_fixed = b_gse_fluctuations[:, 1]
    bz_fixed = b_gse_fluctuations[:, 2]

    # A simple proxy for helicity calculation: Sum(By(t) * Bz(t+1) - Bz(t) * By(t+1))
    # A positive value suggests right-hand polarization, negative suggests left-hand.
    helicity_proxy_fixed = np.sum(by_fixed[:-1] * bz_fixed[1:] - bz_fixed[:-1] * by_fixed[1:])
    
    print(f"\nPerpendicular components assumed to be By and Bz.")
    print(f"Helicity Proxy (Fixed Frame) = {helicity_proxy_fixed:.2f}")
    print("-" * 50)

    # --- Part 2: Using the Field-Aligned Coordinate system (Proper Method) ---
    print("\n--- Part 2: Analysis in Field-Aligned Frame (FAC) ---")
    print("This method finds the true perpendicular plane before calculation.")
    
    # The function prints its step-by-step calculations
    b_par, b_p1, b_p2 = calculate_field_aligned_components(b_gse)

    # Calculate helicity proxy using the correct perpendicular components
    helicity_proxy_fac = np.sum(b_p1[:-1] * b_p2[1:] - b_p2[:-1] * b_p1[1:])

    print(f"\nHelicity Proxy (FAC Frame) = {helicity_proxy_fac:.2f}")
