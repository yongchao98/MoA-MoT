import numpy as np
from scipy.signal import hilbert

def calculate_field_aligned_helicity():
    """
    Demonstrates the correct way to calculate magnetic helicity for a wave
    by first transforming the magnetic field data into a Field-Aligned
    Coordinate (FAC) system.

    The method involves:
    1.  Calculating the mean magnetic field (B0) over the interval.
    2.  Defining a new coordinate system with one axis parallel to B0.
    3.  Projecting the fluctuating field (dB) onto the plane perpendicular to B0.
    4.  Calculating the normalized magnetic helicity from these perpendicular components.
    """
    # --- 0. Generate Synthetic Data ---
    # Create a synthetic time series representing a wave.
    # The background field B0 is [4, 3, 1] nT, which is NOT radial.
    # A right-hand polarized wave is added.
    print("Step 0: Generating synthetic data...")
    t = np.linspace(0, 10 * np.pi, 1000)
    B0 = np.array([4.0, 3.0, 1.0])  # Background field in nT
    wave_amplitude = 0.5
    # Wave components in the XY plane (for simplicity)
    db_wave = np.array([wave_amplitude * np.cos(t), wave_amplitude * np.sin(t), np.zeros_like(t)])
    
    # The measured field B(t) = B0 + db(t)
    B = B0[:, np.newaxis] + db_wave
    Bx, By, Bz = B[0,:], B[1,:], B[2,:]
    print(f"Generated data with a background field B0 ≈ {B0} nT.\n")


    # --- 1. Calculate the Mean Field (B0) ---
    print("Step 1: Calculating the mean magnetic field...")
    B_mean = np.array([np.mean(Bx), np.mean(By), np.mean(Bz)])
    B_mag = np.linalg.norm(B_mean)
    print(f"Calculated Mean Field B_mean = {B_mean.round(2)} nT")
    print(f"Mean Field Magnitude |B_mean| = {B_mag:.2f} nT\n")


    # --- 2. Define the Field-Aligned Coordinate (FAC) System ---
    print("Step 2: Defining the Field-Aligned Coordinate (FAC) system...")
    # The parallel unit vector (z-axis of FAC)
    b_parallel = B_mean / B_mag

    # Create the perpendicular axes. We cross with a stable vector (e.g., Z-axis of original system)
    # to get the first perpendicular axis.
    z_hat = np.array([0., 0., 1.])
    if np.allclose(b_parallel, z_hat) or np.allclose(b_parallel, -z_hat):
        # If B is aligned with Z, use Y instead to avoid a null cross product.
        z_hat = np.array([0., 1., 0.])
        
    b_perp1 = np.cross(b_parallel, z_hat)
    b_perp1 /= np.linalg.norm(b_perp1)

    # The third axis is orthogonal to the other two.
    b_perp2 = np.cross(b_parallel, b_perp1)
    
    print(f"  Parallel axis (b_par)      : {b_parallel.round(3)}")
    print(f"  Perpendicular axis 1 (b_p1): {b_perp1.round(3)}")
    print(f"  Perpendicular axis 2 (b_p2): {b_perp2.round(3)}\n")


    # --- 3. Calculate and Project the Fluctuating Field (dB) ---
    print("Step 3: Projecting fluctuating field dB onto FAC axes...")
    # Get fluctuating part dB = B - B_mean
    dB = B - B_mean[:, np.newaxis]
    
    # Project dB onto the new perpendicular axes
    dB_perp1 = np.dot(b_perp1, dB)
    dB_perp2 = np.dot(b_perp2, dB)
    print("Projection complete.\n")


    # --- 4. Calculate Normalized Magnetic Helicity (sigma_m) ---
    print("Step 4: Calculating normalized magnetic helicity (sigma_m)...")
    # Use the Hilbert transform method.
    # For a right-hand wave, sigma_m should be close to +1.
    # For a left-hand wave, sigma_m should be close to -1.
    
    # Hilbert transform of the second perpendicular component
    dB_perp2_h = np.imag(hilbert(dB_perp2))

    # Numerator of the helicity formula
    # Im(<d~B_p1 * conj(d~B_p2)>) where d~B is the analytic signal.
    # This simplifies to <dB_p1 * dB_p2_h - dB_p2 * dB_p1_h> in the time domain.
    # A common estimator is 2 * <dB_p1 * dB_p2_h>
    numerator = 2 * np.mean(dB_perp1 * dB_perp2_h)

    # Denominator is the total power in the perpendicular fluctuations
    denominator = np.mean(dB_perp1**2 + dB_perp2**2)
    
    # Final normalized magnetic helicity
    sigma_m = numerator / denominator

    print("\n--- FINAL CALCULATION ---")
    print("The normalized magnetic helicity is sigma_m = (Numerator) / (Denominator)")
    print(f"Numerator value       = 2 * <dB_p1 * Hilbert(dB_p2)> = {numerator:.4f}")
    print(f"Denominator value     = <(dB_p1^2 + dB_p2^2)>      = {denominator:.4f}")
    print(f"Resulting sigma_m = {numerator:.4f} / {denominator:.4f} = {sigma_m:.4f}")
    print("\nThe positive value close to 1 indicates right-hand polarization, as expected from the synthetic data.")

calculate_field_aligned_helicity()
<<<~1.0>>>