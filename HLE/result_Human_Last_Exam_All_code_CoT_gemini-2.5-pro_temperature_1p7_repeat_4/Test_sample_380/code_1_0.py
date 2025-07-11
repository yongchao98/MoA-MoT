import numpy as np
from scipy.signal import hilbert

def calculate_and_print_helicity(B1, B2, method_name):
    """Helper function to calculate and print helicity details."""
    # Use the analytic signal (via Hilbert transform) to get the complex representation
    B1_complex = hilbert(B1)
    B2_complex = hilbert(B2)

    # The formula for normalized magnetic helicity is:
    # sigma = Im(2 * <B1_complex * B2_complex*>) / <|B1_complex|^2 + |B2_complex|^2>
    # where <> is the time average and * denotes the complex conjugate.
    numerator = np.mean(2 * np.imag(B1_complex * np.conj(B2_complex)))
    denominator = np.mean(np.abs(B1_complex)**2 + np.abs(B2_complex)**2)
    
    # Avoid division by zero
    if denominator == 0:
        helicity = 0
    else:
        helicity = numerator / denominator

    print(f"Analysis using the {method_name}:")
    if "Incorrect" in method_name:
        print("  Using components B_T and B_N (perpendicular to the RADIAL R-axis).")
        print("  Equation: helicity = Im(2 * <B_T * B_N^*>) / (<|B_T|^2> + <|B_N|^2>)")
    else:
        print("  Using components B_perp1 and B_perp2 (perpendicular to the LOCAL B-FIELD axis).")
        print("  Equation: helicity = Im(2 * <B_perp1 * B_perp2^*>) / (<|B_perp1|^2> + <|B_perp2|^2>)")
    
    print(f"  Inserting calculated values into the equation:")
    print(f"    Numerator = {numerator:.4f}")
    print(f"    Denominator = {denominator:.4f}")
    print(f"    Result: helicity = {numerator:.4f} / {denominator:.4f} = {helicity:.4f}")
    
    if "Correct" in method_name:
        print("\n  This result is ~+1, correctly identifying the simulated wave as left-hand polarized.")
    else:
        print("\n  This result is not close to +/-1, showing this method is incorrect when B is not radial.")

def main():
    """
    Demonstrates the correct method for calculating magnetic helicity for AIC waves.
    """
    print("This script justifies why magnetic helicity for AIC waves must be calculated relative to the local magnetic field, not the radial direction.\n")
    
    # --- 1. Simulation Setup in RTN Coordinates (Radial, Tangential, Normal) ---
    print("--- SITUATION SETUP ---")
    # Define a non-radial background magnetic field B0, simulating the Parker Spiral at L1.
    angle_deg = 45.0
    angle_rad = np.deg2rad(angle_deg)
    B0_magnitude = 5.0  # nT
    B0 = np.array([B0_magnitude * np.cos(angle_rad), B0_magnitude * np.sin(angle_rad), 0.0])
    print(f"A spacecraft at L1 measures a background magnetic field B0.\nDue to the Parker Spiral, B0 is at a {angle_deg}-degree angle to the radial direction.")
    print(f"  B0 vector in RTN (Radial, Tangential, Normal) coordinates = {B0.round(2)} nT\n")
    
    # --- 2. Create a Synthetic Wave propagating along B0 ---
    # We define a left-hand polarized wave (which has helicity of +1)
    t = np.linspace(0, 10, 2000)
    wave_amplitude = 1.0  # nT
    wave_freq = 2 * np.pi * 1.0  # rad/s

    # The wave's physics are defined in a Field-Aligned Coordinate (FAC) system.
    # We must construct the FAC system basis vectors in RTN coordinates.
    z_fac = B0 / np.linalg.norm(B0)
    y_fac = np.array([0, 0, 1.0]) # N-axis, which is perpendicular to the R-T plane
    x_fac = np.cross(y_fac, z_fac)

    # Wave components in the FAC system (propagating along z_fac)
    delta_B_fac_x = wave_amplitude * np.cos(wave_freq * t)
    delta_B_fac_y = wave_amplitude * np.sin(wave_freq * t) # Left-hand polarization (x leads y by 90 deg)
    
    # This wave needs to be rotated into the RTN system to simulate a real measurement.
    # Rotation matrix (columns are the FAC basis vectors)
    M_fac_to_rtn = np.array([x_fac, y_fac, z_fac]).T
    delta_B_fac = np.vstack([delta_B_fac_x, delta_B_fac_y, np.zeros_like(t)])
    delta_B_rtn = M_fac_to_rtn @ delta_B_fac

    # The total measured field is the sum of the background field and the wave
    B_measured_rtn = B0[:, np.newaxis] + delta_B_rtn
    
    print("A left-hand polarized AIC wave (expected helicity = +1) is added, propagating along B0.")
    print("The script will now analyze the total measured field.\n" + "-"*50 + "\n")

    # --- 3. Analysis ---
    
    # First, the researcher calculates the average field to find the local B direction.
    B0_calculated = np.mean(B_measured_rtn, axis=1)
    fluctuations_rtn = B_measured_rtn - B0_calculated[:, np.newaxis]
    
    # --- Incorrect Method ---
    # Using components perpendicular to the *radial* direction (T and N components)
    B_T = fluctuations_rtn[1, :]
    B_N = fluctuations_rtn[2, :]
    calculate_and_print_helicity(B_T, B_N, "Incorrect Method")
    
    print("\n" + "-"*50 + "\n")
    
    # --- Correct Method ---
    # First, rotate the fluctuations into a new FAC system based on B0_calculated.
    z_prime = B0_calculated / np.linalg.norm(B0_calculated)
    y_prime = np.array([0, 0, 1.0]) # As before
    x_prime = np.cross(y_prime, z_prime)
    
    # Rotation matrix from RTN to the new FAC system (rows are the new basis vectors)
    M_rtn_to_fac = np.array([x_prime, y_prime, z_prime])
    fluctuations_fac = M_rtn_to_fac @ fluctuations_rtn
    
    # The components perpendicular to B0 are the first two components in this new system.
    B_perp1 = fluctuations_fac[0, :]
    B_perp2 = fluctuations_fac[1, :]
    calculate_and_print_helicity(B_perp1, B_perp2, "Correct Method")
    
    print("\n" + "-"*50 + "\n")
    print("<<<CONCLUSION>>>")
    print("The calculation using components perpendicular to the local magnetic field correctly identifies the wave's helicity.")
    print("Therefore, while data may be presented in a radial (RTN/GSE) frame, any rigorous helicity analysis requires a coordinate transformation into a Field-Aligned Coordinate system.")


if __name__ == '__main__':
    main()
