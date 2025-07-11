import numpy as np
import matplotlib.pyplot as plt

def calculate_s11_db(R, L, C, freq_hz, Z0=50):
    """
    Calculates the S11 (reflection coefficient) in dB for a series RLC circuit.
    This simulates what a VNA measures.

    Args:
        R (float): Resistance in Ohms.
        L (float): Inductance in Henries.
        C (float): Capacitance in Farads.
        freq_hz (np.ndarray): Array of frequencies in Hz.
        Z0 (int): Characteristic impedance of the measurement system (usually 50 Ohms).

    Returns:
        np.ndarray: S11 value in dB for each frequency.
    """
    omega = 2 * np.pi * freq_hz
    # Impedance of the RLC circuit
    Z_l = 1j * omega * L
    Z_c = 1 / (1j * omega * C)
    Z_coil = R + Z_l + Z_c
    
    # Calculate reflection coefficient (Gamma)
    gamma = (Z_coil - Z0) / (Z_coil + Z0)
    
    # Calculate S11 in dB
    s11_db = 20 * np.log10(np.abs(gamma))
    return s11_db

# --- Simulation Parameters ---
# We will design a coil to be resonant around 64 MHz (for a 1.5T MRI)
# f = 1 / (2*pi*sqrt(L*C))
L_coil = 300e-9  # Inductance = 300 nH
C_coil = 20.7e-12 # Capacitance = 20.7 pF
Z0 = 50           # System Impedance = 50 Ohms

# --- Case 1: Normal, Resonating Coil (High-Q) ---
# A high-Q coil has very low resistance.
R_high_Q = 0.5    # Resistance = 0.5 Ohms

# --- Case 2: Actively Decoupled Coil (Low-Q) ---
# When decoupled, a PIN diode might add significant resistance to the circuit,
# "spoiling" the quality factor (Q) and killing the resonance.
R_decoupled = 500 # Resistance = 500 Ohms

# Frequency range for our simulation
frequencies = np.linspace(50e6, 80e6, 501) # 50 MHz to 80 MHz

# --- Calculate S11 for both cases ---
s11_high_Q = calculate_s11_db(R_high_Q, L_coil, C_coil, frequencies, Z0)
s11_decoupled = calculate_s11_db(R_decoupled, L_coil, C_coil, frequencies, Z0)
resonant_freq_mhz = 1 / (2 * np.pi * np.sqrt(L_coil * C_coil)) / 1e6

print("--- Simulation Parameters ---")
print(f"Coil Inductance (L): {L_coil * 1e9:.1f} nH")
print(f"Coil Capacitance (C): {C_coil * 1e12:.2f} pF")
print(f"Calculated Resonant Frequency: {resonant_freq_mhz:.2f} MHz\n")

print(f"Case 1 (High-Q) Resistance (R): {R_high_Q} Ohms")
print("This simulates the coil in its 'Receive' state, showing a clear resonance dip.\n")

print(f"Case 2 (Decoupled) Resistance (R): {R_decoupled} Ohms")
print("This simulates the coil in its default 'Detuned' state. The resonance is completely flattened.\n")

# --- Plot the results ---
plt.figure(figsize=(10, 6))
plt.plot(frequencies / 1e6, s11_high_Q, label=f'Resonating Coil (R={R_high_Q}Ω)', linewidth=2)
plt.plot(frequencies / 1e6, s11_decoupled, label=f'Decoupled Coil (R={R_decoupled}Ω)', linewidth=2, linestyle='--')
plt.title('VNA Measurement (S11) of an MRI Coil')
plt.xlabel('Frequency (MHz)')
plt.ylabel('S11 - Reflection Coefficient (dB)')
plt.grid(True)
plt.legend()
plt.ylim(-40, 1) # Set y-axis limits for better visualization
plt.axvline(x=resonant_freq_mhz, color='r', linestyle=':', label=f'f_res = {resonant_freq_mhz:.2f} MHz')
plt.legend()
plt.show()
