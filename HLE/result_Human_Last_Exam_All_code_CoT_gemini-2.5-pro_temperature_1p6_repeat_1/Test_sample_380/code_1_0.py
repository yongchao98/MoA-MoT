import numpy as np
from scipy.signal import hilbert

def calculate_helicity():
    """
    This function demonstrates the correct way to calculate magnetic helicity
    by comparing an approximate method with the physically correct Field-Aligned
    Coordinate (FAC) method.
    """
    # 1. SIMULATE SPACECRAFT DATA
    # We will create a synthetic magnetic field as a spacecraft at L1 might see it.
    # The coordinate system is like GSE, where X is roughly Sun-Earth (radial).
    
    # Define a mean magnetic field (B0) that follows a Parker Spiral (not purely radial).
    # Here, Bx and By components are equal, representing a 45-degree spiral angle.
    B0_nT = np.array([3.5, 3.5, 0.5])  # In nanotesla (nT)

    # Define the properties of a pure Alfvén Ion Cyclotron (AIC) wave.
    # We will create a purely left-hand polarized wave.
    # The wave exists in the plane perpendicular to the mean field B0.
    
    # First, define the Field-Aligned Coordinate (FAC) system based on B0.
    B0_mag = np.linalg.norm(B0_nT)
    z_fac = B0_nT / B0_mag  # The axis parallel to B0
    
    # Create a perpendicular axis by crossing with a stable direction (e.g., Z-axis)
    # This defines the perpendicular plane.
    y_fac_temp = np.cross(z_fac, [0, 0, 1])
    y_fac = y_fac_temp / np.linalg.norm(y_fac_temp)
    x_fac = np.cross(y_fac, z_fac) # This completes the right-handed system {x_fac, y_fac, z_fac}

    # Generate the synthetic wave signal in the FAC perpendicular plane {x_fac, y_fac}
    t = np.linspace(0, 50, 2000)
    amplitude_nT = 0.5
    frequency_hz = 0.3
    # B_y leads B_x by 90 deg for a left-hand polarized wave (positive helicity).
    # Let's use the standard helicity formula where B1*H(B2) - B2*H(B1) gives negative
    # for right-hand and positive for left-hand pol.
    # B1 = amp * cos(wt), B2 = amp * sin(wt) should give left-hand.
    db_perp1 = amplitude_nT * np.cos(2 * np.pi * frequency_hz * t)
    db_perp2 = amplitude_nT * np.sin(2 * np.pi * frequency_hz * t)

    # Transform the wave from FAC back to the original (GSE) coordinates.
    # This is what the spacecraft instrument would measure.
    # The total fluctuating field db = db_perp1 * x_fac_vec + db_perp2 * y_fac_vec
    db_gse = np.outer(db_perp1, x_fac) + np.outer(db_perp2, y_fac)
    
    # The total measured magnetic field is the sum of the mean field and the wave.
    B_total_gse = B0_nT + db_gse
    Bx, By, Bz = B_total_gse[:, 0], B_total_gse[:, 1], B_total_gse[:, 2]

    # --- START ANALYSIS ---
    # Now, we pretend we are a scientist who just has the time series Bx, By, Bz.

    # 2. APPROXIMATE METHOD: USING Y AND Z COMPONENTS
    # This method assumes the radial direction is X and uses Y,Z as the perpendicular plane.
    # As you suspected, this is physically inaccurate because B0 is not along X.
    print("--- Method 1: Approximation using Y and Z components ---")
    db_y_approx = By - np.mean(By)
    db_z_approx = Bz - np.mean(Bz)
    
    # Use the Hilbert transform to calculate the normalized magnetic helicity.
    db_y_hilbert = np.imag(hilbert(db_y_approx))
    db_z_hilbert = np.imag(hilbert(db_z_approx))

    num_approx = db_y_approx * db_z_hilbert - db_z_approx * db_y_hilbert
    den_approx = db_y_approx**2 + db_z_approx**2
    # To avoid division by zero, replace any zeros in the denominator.
    den_approx[den_approx == 0] = 1e-9
    
    # The final equation is the mean of the numerator divided by the mean of the denominator.
    final_num_approx = np.mean(num_approx)
    final_den_approx = np.mean(den_approx)
    helicity_approx = final_num_approx / final_den_approx

    print(f"Components used for 'perpendicular' plane: Fluctuating By, Bz")
    print(f"Final equation: {final_num_approx:.4f} / {final_den_approx:.4f}")
    print(f"Calculated Normalized Helicity = {helicity_approx:.4f}")
    print("(This result is small because it incorrectly mixes parallel and perpendicular components.)\n")


    # 3. CORRECT METHOD: USING A FIELD-ALIGNED COORDINATE (FAC) SYSTEM
    # This is the physically correct way to perform the analysis.
    print("--- Method 2: Correct analysis using Field-Aligned Coordinates ---")
    
    # First, calculate the mean magnetic field from the data.
    B_mean_calc = np.mean(B_total_gse, axis=0)
    
    # Define the FAC unit vectors based on this calculated mean field.
    z_hat = B_mean_calc / np.linalg.norm(B_mean_calc)
    y_hat_temp = np.cross(z_hat, [0, 0, 1])
    y_hat = y_hat_temp / np.linalg.norm(y_hat_temp)
    x_hat = np.cross(y_hat, z_hat)

    # Project the fluctuating magnetic field onto the new FAC basis.
    db_gse_calc = B_total_gse - B_mean_calc
    db_correct_perp1 = np.dot(db_gse_calc, x_hat)
    db_correct_perp2 = np.dot(db_gse_calc, y_hat)
    
    # Calculate helicity using the two correct perpendicular components.
    db_perp1_hilbert = np.imag(hilbert(db_correct_perp1))
    db_perp2_hilbert = np.imag(hilbert(db_correct_perp2))

    num_correct = db_correct_perp1 * db_perp2_hilbert - db_correct_perp2 * db_perp1_hilbert
    den_correct = db_correct_perp1**2 + db_correct_perp2**2
    den_correct[den_correct == 0] = 1e-9

    final_num_correct = np.mean(num_correct)
    final_den_correct = np.mean(den_correct)
    helicity_correct = final_num_correct / final_den_correct
    
    print(f"Mean Field B0 (from data) = {B_mean_calc[0]:.2f}, {B_mean_calc[1]:.2f}, {B_mean_calc[2]:.2f} nT")
    print(f"Components used for perpendicular plane: Projections onto vectors perpendicular to B0.")
    print(f"Final equation: {final_num_correct:.4f} / {final_den_correct:.4f}")
    print(f"Calculated Normalized Helicity = {helicity_correct:.4f}")
    print("(Note: A value near +1.0 indicates pure left-hand polarization, as expected from our input wave.)\n")


if __name__ == '__main__':
    calculate_helicity()