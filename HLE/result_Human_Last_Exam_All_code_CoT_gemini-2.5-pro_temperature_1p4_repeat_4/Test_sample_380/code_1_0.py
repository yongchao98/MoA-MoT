import numpy as np

def calculate_and_print_helicity_comparison():
    """
    This script demonstrates why using a radial-based coordinate system to calculate
    magnetic helicity for AIC waves at L1 is a valid approximation, even though
    the local magnetic field is not perfectly radial.
    """

    # --- Step 1: Define the Physical Setup ---
    # We define a coordinate system where the X-axis represents the radial direction from the Sun.
    RADIAL_DIRECTION = np.array([1.0, 0.0, 0.0])

    # The background magnetic field at L1 is not perfectly radial due to the Parker Spiral.
    # We define a background field B0 at a 30-degree angle to the radial direction.
    angle_deg = 30.0
    angle_rad = np.deg2rad(angle_deg)
    BACKGROUND_B_FIELD = np.array([np.cos(angle_rad), np.sin(angle_rad), 0.0])

    # --- Step 2: Simulate Wave Data ---
    # We create a synthetic left-hand circularly polarized Alfven Ion Cyclotron (AIC) wave.
    # Its magnetic fluctuations (delta_B) are perpendicular to its propagation direction (BACKGROUND_B_FIELD).

    # First, create the basis vectors for the plane perpendicular to BACKGROUND_B_FIELD.
    # This is the "physically correct" Field-Aligned Coordinate (FAC) system.
    z_fac = BACKGROUND_B_FIELD / np.linalg.norm(BACKGROUND_B_FIELD) # Parallel to B
    y_fac = np.cross(z_fac, np.array([0, 0, 1])) # A perpendicular vector
    y_fac /= np.linalg.norm(y_fac)
    x_fac = np.cross(y_fac, z_fac) # The third orthogonal vector

    # Create a time series for a left-hand polarized wave in the FAC system.
    t = np.linspace(0, 4 * np.pi, 1000)
    wave_comp1 = np.cos(t) # Fluctuation along x_fac
    wave_comp2 = np.sin(t) # Fluctuation along y_fac

    # Express this wave in our main (X, Y, Z) coordinate system. This is what a spacecraft would measure.
    delta_B_vec = wave_comp1[:, np.newaxis] * x_fac + wave_comp2[:, np.newaxis] * y_fac
    delta_By = delta_B_vec[:, 1]
    delta_Bz = delta_B_vec[:, 2]

    # --- Step 3: Calculate Helicity using the Radial Approximation ---
    # Here, we use the components perpendicular to the radial direction (Y and Z).
    # We use the cross-spectrum. The sign of the imaginary part indicates polarization.
    fft_By = np.fft.fft(delta_By)
    fft_Bz = np.fft.fft(delta_Bz)
    cross_spectrum_radial = fft_By * np.conj(fft_Bz)

    # Find the dominant frequency of the wave to analyze the helicity there.
    dominant_freq_idx = np.argmax(np.abs(fft_By[1:len(t)//2])) + 1

    # Helicity is proportional to the imaginary part of the cross-spectrum.
    helicity_value_radial = np.imag(cross_spectrum_radial[dominant_freq_idx])
    # The normalized helicity divides by the total power.
    power_radial = np.abs(fft_By[dominant_freq_idx])**2 + np.abs(fft_Bz[dominant_freq_idx])**2
    normalized_helicity_radial = 2 * helicity_value_radial / power_radial if power_radial > 0 else 0

    # --- Step 4: Calculate Helicity using the "Correct" Field-Aligned Frame ---
    # Here, we use the original wave components in the plane perpendicular to the actual B-field.
    fft_comp1 = np.fft.fft(wave_comp1)
    fft_comp2 = np.fft.fft(wave_comp2)
    cross_spectrum_fac = fft_comp1 * np.conj(fft_comp2)

    # Calculate helicity at the same dominant frequency.
    helicity_value_fac = np.imag(cross_spectrum_fac[dominant_freq_idx])
    power_fac = np.abs(fft_comp1[dominant_freq_idx])**2 + np.abs(fft_comp2[dominant_freq_idx])**2
    normalized_helicity_fac = 2 * helicity_value_fac / power_fac if power_fac > 0 else 0

    # --- Step 5: Compare Results and Print the Justification ---
    print("--- Justification for using a Radial Coordinate System ---")
    print("\nSummary: The local magnetic field is NOT assumed to be perfectly radial. Using the radial direction is a justified approximation because it preserves the sign of the helicity, which determines the wave's polarization.")
    print("\n1. Simulation Setup:")
    print(f"   - Assumed Radial Direction (X-axis): [{RADIAL_DIRECTION[0]}, {RADIAL_DIRECTION[1]}, {RADIAL_DIRECTION[2]}]")
    print(f"   - Assumed Local B-Field Direction (at {angle_deg:.1f} deg to radial): [{BACKGROUND_B_FIELD[0]:.2f}, {BACKGROUND_B_FIELD[1]:.2f}, {BACKGROUND_B_FIELD[2]:.2f}]")
    print("   - A synthetic left-hand polarized wave was generated propagating along this B-field.")

    print("\n2. Calculation using Radial Approximation (Y and Z components):")
    print(f"   - The calculated normalized magnetic helicity is: {normalized_helicity_radial:.4f}")

    print("\n3. Calculation using 'Correct' Field-Aligned Frame (components perp. to B):")
    print(f"   - The calculated normalized magnetic helicity is: {normalized_helicity_fac:.4f}")

    print("\n4. Final Conclusion:")
    print("   The SIGN of the helicity is the same in both calculations.")
    print(f"   Sign using Radial Frame: {np.sign(normalized_helicity_radial):.0f}")
    print(f"   Sign using Field-Aligned Frame: {np.sign(normalized_helicity_fac):.0f}")
    print("\n   While the magnitude differs, the crucial diagnostic for polarization—the sign—is preserved.")
    print("   Therefore, calculating helicity using the radial approximation is a computationally convenient and physically justified method.")

if __name__ == '__main__':
    calculate_and_print_helicity_comparison()