import numpy as np
import matplotlib.pyplot as plt

def model_s11_dip(freq_array, f0, Q, coupling_factor):
    """
    Models the magnitude of the S11 reflection coefficient for a resonator.
    A deeper dip corresponds to better coupling and higher Q.
    
    Args:
        freq_array (np.array): Array of frequencies to calculate S11 for.
        f0 (float): The center resonance frequency.
        Q (float): The Quality Factor of the resonator.
        coupling_factor (float): A value from 0 to 1 representing coupling depth.

    Returns:
        np.array: The magnitude of S11 (from 0 to 1).
    """
    # A Lorentzian function is a good model for a resonance shape.
    # The bandwidth (FWHM) of the resonance is f0 / Q.
    lorentzian = 1 / (1 + ( (freq_array - f0) / (f0 / (2 * Q)) )**2)
    
    # S11 magnitude is 1 (total reflection) minus the dip caused by the resonance.
    s11_mag = 1.0 - coupling_factor * lorentzian
    return s11_mag

def main():
    # --- Simulation Parameters ---
    # Define a frequency range around a typical 1.5T MRI frequency (proton Larmor freq)
    f_start_mhz = 60
    f_stop_mhz = 68
    points = 801
    freqs_mhz = np.linspace(f_start_mhz, f_stop_mhz, points)
    
    # Coil's intrinsic properties
    f0_mhz = 63.8  # Target resonance frequency in MHz
    coupling = 0.8 # Assume good coupling from the probe, determines max dip depth
    
    # --- Scenario 1: Tuned Coil (High-Q, as in the scanner) ---
    # A high Q-factor means a sharp, narrow resonance.
    Q_tuned = 120
    s11_tuned_mag = model_s11_dip(freqs_mhz, f0_mhz, Q_tuned, coupling)
    s11_tuned_db = 20 * np.log10(s11_tuned_mag)
    
    # --- Scenario 2: Detuned Coil (Low-Q, default state on the bench) ---
    # The detuning circuit adds resistance, killing the Q-factor.
    Q_detuned = 5
    s11_detuned_mag = model_s11_dip(freqs_mhz, f0_mhz, Q_detuned, coupling)
    s11_detuned_db = 20 * np.log10(s11_detuned_mag)

    # --- Plotting the results ---
    plt.figure(figsize=(12, 7))
    plt.plot(freqs_mhz, s11_tuned_db, label=f'Tuned Coil (Q={Q_tuned}) - Resonance Visible', linewidth=2)
    plt.plot(freqs_mhz, s11_detuned_db, label=f'Detuned Coil (Q={Q_detuned}) - Resonance "Invisible"', linewidth=2, linestyle='--')
    
    plt.title('Simulated VNA Measurement (S11) of an MRI Coil')
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('S11 (dB)')
    plt.grid(True, which='both', linestyle=':')
    plt.ylim(-20, 1)
    plt.legend()
    
    print("This plot simulates the VNA measurement for a tuned vs. detuned MRI coil.")
    print("The 'Tuned' state shows a sharp, deep resonance dip, as expected during a scan.")
    print("The 'Detuned' state, caused by safety circuitry, has such a low Q-factor")
    print("that the dip is broad and shallow, making it practically invisible.\n")
    
    # --- Output the equation parameters for the visible resonance ---
    print("------------------------------------------------------------------")
    print("Model Parameters for the Visible Resonance (Tuned Coil) Curve:")
    print("The curve is based on a Lorentzian dip function model.")
    
    # Bandwidth (Full Width at Half Maximum) of the resonance peak
    bandwidth_mhz = f0_mhz / Q_tuned
    
    print(f"Center Resonance Frequency (f0): {f0_mhz} MHz")
    print(f"Quality Factor (Q): {Q_tuned}")
    print(f"Resulting Bandwidth (f0/Q): {bandwidth_mhz:.4f} MHz")
    print(f"Coupling Factor (determines max dip depth): {coupling}")
    print("------------------------------------------------------------------")
    
    plt.show()

if __name__ == '__main__':
    main()