import numpy as np

def calculate_s11_db(R, L, C, freq_hz, z0=50):
    """
    Calculates the S11 reflection coefficient in dB for a series RLC circuit.
    This simulates what a Vector Network Analyzer (VNA) measures.

    Args:
        R (float): Resistance in Ohms.
        L (float): Inductance in Henries.
        C (float): Capacitance in Farads.
        freq_hz (np.ndarray): Array of frequencies in Hertz.
        z0 (int): Characteristic impedance of the system (usually 50 Ohms).

    Returns:
        np.ndarray: S11 values in dB for each frequency.
    """
    omega = 2 * np.pi * freq_hz
    # Impedance of the series RLC circuit (the coil)
    z_coil = R + 1j * (omega * L - 1 / (omega * C))
    # Reflection coefficient Gamma
    gamma = (z_coil - z0) / (z_coil + z0)
    # S11 in dB
    s11_db = 20 * np.log10(np.abs(gamma))
    return s11_db

# --- Simulation Parameters ---
# Let's assume a coil designed for a 3T MRI scanner (proton frequency ~127.7 MHz)
L_coil = 100e-9  # Inductance = 100 nH
C_coil = 15.53e-12 # Capacitance = 15.53 pF (chosen to resonate near 127.7 MHz)
f0 = 1 / (2 * np.pi * np.sqrt(L_coil * C_coil)) / 1e6 # Calculate resonant freq in MHz

# Define frequency range for the simulation sweep
frequencies = np.linspace(120e6, 135e6, 1001) # 120 MHz to 135 MHz

# --- Scenario 1: High-Q Coil (unloaded, isolated on the bench) ---
# A high-quality coil has very low intrinsic resistance.
R_high_q = 0.5  # Ohms (low resistance -> high Q)

s11_high_q = calculate_s11_db(R_high_q, L_coil, C_coil, frequencies)
min_s11_high_q = np.min(s11_high_q)
res_freq_high_q = frequencies[np.argmin(s11_high_q)] / 1e6

print("--- Scenario 1: High-Q (Undamped) Coil ---")
print("This simulates an isolated coil with a visible resonance.")
print(f"Coil parameters: R = {R_high_q} Ohms, L = {L_coil*1e9:.2f} nH, C = {C_coil*1e12:.2f} pF")
print(f"The theoretical resonant frequency is: {f0:.2f} MHz")
print(f"The simulated resonance is a sharp dip of {min_s11_high_q:.2f} dB at {res_freq_high_q:.2f} MHz.\n")


# --- Scenario 2: Low-Q Coil (damped by a preamplifier) ---
# The low input impedance of the preamplifier effectively adds a lot of resistance
# to the resonant circuit, "damping" the resonance.
R_low_q = 45.0  # Ohms (high effective resistance -> low Q)

s11_low_q = calculate_s11_db(R_low_q, L_coil, C_coil, frequencies)
min_s11_low_q = np.min(s11_low_q)
res_freq_low_q = frequencies[np.argmin(s11_low_q)] / 1e6

print("--- Scenario 2: Low-Q (Damped) Coil ---")
print("This simulates the same coil when damped by its preamplifier.")
print(f"Coil parameters: R = {R_low_q} Ohms, L = {L_coil*1e9:.2f} nH, C = {C_coil*1e12:.2f} pF")
print("Note: The 'R' here is the *effective* series resistance, representing the damping.")
print(f"The simulated resonance is now a very shallow dip of {min_s11_low_q:.2f} dB at {res_freq_low_q:.2f} MHz.")
print("\nThis extremely shallow dip would be very difficult or impossible to see on a VNA,")
print("making the coil appear non-resonant, even though it is fully functional.")
