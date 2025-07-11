import numpy as np

def simulate_s11_dip(freq_mhz, f0, Q, coupling_factor):
    """
    Simulates the S11 reflection coefficient (in dB) for a resonator.
    A simple Lorentzian model is used for the absorbed power.
    """
    # The term (f0 / (2 * Q)) is the half-width at half-maximum of the resonance power peak.
    normalized_freq_sq = ((freq_mhz - f0) / (f0 / (2 * Q)))**2
    
    # Reflected power is modeled as 1 minus the absorbed power (Lorentzian dip)
    reflected_power_ratio = 1.0 - coupling_factor / (1.0 + normalized_freq_sq)
    
    # S11 in dB is 10 * log10(Power Ratio)
    s11_db = 10 * np.log10(reflected_power_ratio)
    
    # Clamp values to a realistic VNA floor to avoid -inf from log10(0)
    s11_db = np.maximum(s11_db, -40.0)
    
    return s11_db

# --- Simulation Parameters ---
f0 = 63.87  # Resonant frequency in MHz (for a 1.5T MRI)
Q_high = 250  # Q-factor of a high-quality, unloaded resonant coil
Q_low = 30    # Q-factor of a heavily damped or poorly-coupled coil
# This represents how much the probe couples to the coil. A value of 1.0
# would mean perfect energy transfer at resonance. We'll assume moderate coupling.
coupling_factor = 0.8 

# Frequency range for our simulated VNA sweep
frequencies = np.linspace(61, 67, 61) # Sweep from 61 to 67 MHz in 0.1 MHz steps

# --- Run Simulation ---
s11_high_q = simulate_s11_dip(frequencies, f0, Q_high, coupling_factor)
s11_low_q = simulate_s11_dip(frequencies, f0, Q_low, coupling_factor)

# --- Print Results ---
print("--- S11 Simulation: High-Q vs. Low-Q Coil ---")
print(f"Parameters used:")
print(f"  Resonant Frequency (f0): {f0} MHz")
print(f"  High Q Value: {Q_high}")
print(f"  Low Q Value: {Q_low}")
print("-" * 50)
print(f"{'Frequency (MHz)':<20} {'High-Q S11 (dB)':<20} {'Low-Q S11 (dB)':<20}")
print("-" * 50)

for i in range(len(frequencies)):
    # Highlight the resonant frequency
    marker = "<- Resonance" if abs(frequencies[i] - f0) < 0.05 else ""
    print(f"{frequencies[i]:<20.2f} {s11_high_q[i]:<20.3f} {s11_low_q[i]:<20.3f} {marker}")

print("\n--- Observation ---")
print("Notice the High-Q coil has a very deep, sharp dip near resonance (~ -7.0 dB).")
print("The Low-Q coil has a very shallow, broad dip (~ -0.9 dB), which could be hard to see on a real VNA.")
