import numpy as np
from scipy.signal import hilbert

def calculate_helicity():
    """
    This function demonstrates the rigorous method for calculating magnetic helicity
    for a wave in a non-radial background magnetic field.
    """
    # --- 1. Generate Synthetic Spacecraft Data ---
    # We will create a synthetic magnetic field time series as if measured by a
    # spacecraft at L1 in a standard coordinate system (e.g., GSE, where X is Sun-Earth).

    # Define the properties of our synthetic solar wind
    # The average magnetic field (B0) is non-radial, consistent with the Parker Spiral.
    # B0 = [Bx, By, Bz] in nT. Here, Bx is the radial component.
    B0_vec = np.array([3.0, 3.0, 0.5])
    
    # Define the properties of a pure AIC wave propagating along B0
    # This wave is transverse, so its fluctuations are perpendicular to B0.
    # We will create a left-hand circularly polarized wave.
    num_points = 2048
    fs = 10.0  # Sampling frequency in Hz
    t = np.arange(num_points) / fs
    wave_freq = 0.5  # Wave frequency in Hz
    wave_amp = 0.5   # Wave amplitude in nT

    # To create the wave correctly, we first define a basis aligned with B0.
    e_parallel = B0_vec / np.linalg.norm(B0_vec)
    # Create a reference vector (e.g., Z-axis) to define the perpendicular plane.
    ref_vec = np.array([0.0, 0.0, 1.0])
    e_perp1 = np.cross(e_parallel, ref_vec)
    e_perp1 /= np.linalg.norm(e_perp1)
    e_perp2 = np.cross(e_parallel, e_perp1)

    # Create the wave fluctuations in the field-aligned coordinate (FAC) system.
    # For a left-hand polarized wave, the components are cos(wt) and sin(wt).
    deltaB_perp1_fac = wave_amp * np.cos(2 * np.pi * wave_freq * t)
    deltaB_perp2_fac = wave_amp * np.sin(2 * np.pi * wave_freq * t)
    
    # The fluctuation vector at each time step in the FAC system
    # deltaB_fac(t) = deltaB_perp1_fac(t) * e_perp1 + deltaB_perp2_fac(t) * e_perp2
    # We transform this vector back to the original GSE system for our "measurement".
    deltaB_gse = np.outer(deltaB_perp1_fac, e_perp1) + np.outer(deltaB_perp2_fac, e_perp2)

    # The final "measured" magnetic field is the sum of the mean field and the wave.
    B_gse = B0_vec + deltaB_gse
    
    print("--- Step 1: Synthetic Data Generation ---")
    print(f"Generated synthetic data with a non-radial average magnetic field.")
    print(f"Average Field B0 (GSE) = {B0_vec[0]:.2f} X, {B0_vec[1]:.2f} Y, {B0_vec[2]:.2f} Z (nT)\n")

    # --- 2. Analyze the "Measured" Data ---
    # Now, we pretend we only have the B_gse data and want to find the helicity.
    
    # First, calculate the mean field from the time series.
    B0_mean_calculated = np.mean(B_gse, axis=0)
    
    # Calculate the magnetic field fluctuations.
    deltaB = B_gse - B0_mean_calculated

    print("--- Step 2: Field-Aligned Coordinate Transformation ---")
    print(f"Calculated average field from data: {B0_mean_calculated[0]:.2f} X, {B0_mean_calculated[1]:.2f} Y, {B0_mean_calculated[2]:.2f} Z (nT)")
    print("This field is not purely radial (not just in X), so we must transform coordinates.")

    # Define the new field-aligned coordinate (FAC) system based on the calculated mean field.
    e_parallel_calc = B0_mean_calculated / np.linalg.norm(B0_mean_calculated)
    e_perp1_calc = np.cross(e_parallel_calc, ref_vec)
    e_perp1_calc /= np.linalg.norm(e_perp1_calc)
    e_perp2_calc = np.cross(e_parallel_calc, e_perp1_calc)
    
    print("\nNew Field-Aligned Basis Vectors:")
    print(f"e_parallel = [{e_parallel_calc[0]:.3f}, {e_parallel_calc[1]:.3f}, {e_parallel_calc[2]:.3f}]")
    print(f"e_perp1    = [{e_perp1_calc[0]:.3f}, {e_perp1_calc[1]:.3f}, {e_perp1_calc[2]:.3f}]")
    print(f"e_perp2    = [{e_perp2_calc[0]:.3f}, {e_perp2_calc[1]:.3f}, {e_perp2_calc[2]:.3f}]\n")

    # Project the magnetic field fluctuations onto the new perpendicular basis vectors.
    # These are the components that should be used for the helicity calculation.
    deltaB_perp1 = np.dot(deltaB, e_perp1_calc)
    deltaB_perp2 = np.dot(deltaB, e_perp2_calc)

    # --- 3. Calculate Normalized Magnetic Helicity ---
    # A common definition of normalized magnetic helicity (or circular polarization)
    # uses the Hilbert transform.
    # sigma_m = 2 * <dB_perp1 * H(dB_perp2)> / <(dB_perp1^2 + dB_perp2^2)>
    # where H() is the Hilbert transform and <...> denotes a time average.
    # The sign indicates polarization: +1 for left-hand, -1 for right-hand.

    # Use scipy's hilbert function. We need the imaginary part.
    deltaB_perp2_h = hilbert(deltaB_perp2).imag
    
    # Calculate the numerator and denominator for the helicity formula
    numerator = 2 * np.mean(deltaB_perp1 * deltaB_perp2_h)
    denominator = np.mean(deltaB_perp1**2 + deltaB_perp2**2)
    
    # Avoid division by zero if there are no fluctuations
    if denominator == 0:
        sigma_m = 0.0
    else:
        sigma_m = numerator / denominator

    print("--- Step 3: Helicity Calculation ---")
    print("Calculated the two perpendicular fluctuation components in the new system.")
    print(f"Calculated Normalized Magnetic Helicity (sigma_m) = {sigma_m:.4f}")
    print("\n(Note: A value near +1 indicates left-hand polarization, as generated.)")
    
    return sigma_m

if __name__ == '__main__':
    helicity_value = calculate_helicity()
    # The final answer format requested by the user.
    print(f"\n<<<Final Answer>>>")
    print(f"<<<{helicity_value:.4f}>>>")
