import numpy as np
from scipy import signal

def calculate_magnetic_helicity_demo():
    """
    This function demonstrates the calculation of normalized magnetic helicity
    for a synthetic left-hand polarized wave, similar to an AIC wave.
    """
    # 1. Define wave and signal parameters
    fs = 100.0  # Sampling frequency in Hz
    duration = 10  # Signal duration in seconds
    n_points = int(fs * duration) # Total number of data points
    t = np.linspace(0, duration, n_points, endpoint=False) # Time vector

    wave_freq = 5.0  # Frequency of the AIC wave in Hz
    amplitude = 1.0  # Amplitude of the wave
    noise_power = 0.05 # Add some noise for realism

    # 2. Generate synthetic data for a left-hand (LH) polarized wave
    # For an LH wave propagating along X, By leads Bz by 90 degrees.
    # Bz(t) = A * sin(wt), By(t) = A * cos(wt)
    # We add some Gaussian noise to both components.
    b_y = amplitude * np.cos(2 * np.pi * wave_freq * t) + np.random.normal(scale=np.sqrt(noise_power), size=t.shape)
    b_z = amplitude * np.sin(2 * np.pi * wave_freq * t) + np.random.normal(scale=np.sqrt(noise_power), size=t.shape)

    # 3. Calculate power spectra and cross-spectrum
    # We use Welch's method for a robust spectral estimate.
    nperseg = 256 # Window size for Welch's method
    
    # Power Spectral Density of By (P_yy)
    freqs, p_yy = signal.welch(b_y, fs, nperseg=nperseg)
    # Power Spectral Density of Bz (P_zz)
    _, p_zz = signal.welch(b_z, fs, nperseg=nperseg)
    # Cross Spectral Density between By and Bz (P_yz)
    _, p_yz = signal.csd(b_y, b_z, fs, nperseg=nperseg)

    # 4. Find the spectral values at the wave's frequency
    # Find the index in the frequency array closest to our wave frequency
    wave_idx = np.argmin(np.abs(freqs - wave_freq))

    # Extract the spectral values at that frequency
    p_yy_wave = p_yy[wave_idx]
    p_zz_wave = p_zz[wave_idx]
    p_yz_wave = p_yz[wave_idx] # This is a complex number

    # 5. Calculate the normalized magnetic helicity
    # The equation is: sigma = 2 * Im(P_yz) / (P_yy + P_zz)
    # The imaginary part of the CSD gives the phase relationship.
    imag_p_yz = np.imag(p_yz_wave)
    sum_p_yy_p_zz = p_yy_wave + p_zz_wave
    
    # Avoid division by zero if the signal is pure noise with no power
    if sum_p_yy_p_zz > 0:
        sigma = (2 * imag_p_yz) / sum_p_yy_p_zz
    else:
        sigma = 0

    # Print the results and the components of the final equation
    print("--- Magnetic Helicity Calculation Demo ---")
    print(f"Analysis performed at frequency: {freqs[wave_idx]:.2f} Hz\n")
    print("Equation: sigma = 2 * Im(P_yz) / (P_yy + P_zz)\n")
    print("Component values:")
    print(f"  P_yy (Power in Y component) = {p_yy_wave:.4f}")
    print(f"  P_zz (Power in Z component) = {p_zz_wave:.4f}")
    print(f"  P_yz (Cross-power Y-Z)     = {p_yz_wave:.4f}")
    print(f"  Im(P_yz)                     = {imag_p_yz:.4f}\n")
    print("Final Calculation:")
    print(f"  sigma = (2 * {imag_p_yz:.4f}) / ({p_yy_wave:.4f} + {p_zz_wave:.4f})")
    print(f"  sigma = {2 * imag_p_yz:.4f} / {sum_p_yy_p_zz:.4f}")
    print(f"  Normalized Magnetic Helicity (sigma) = {sigma:.4f}\n")
    print("Note: A value close to +1 indicates left-hand polarization, as expected.")

# Execute the demonstration
calculate_magnetic_helicity_demo()
<<<+1>>>