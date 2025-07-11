import numpy as np

def calculate_s11_db(R, L, C, Z0, frequency):
    """
    Calculates the S11 (reflection coefficient) in dB for a series RLC circuit.

    Args:
        R (float): Resistance in Ohms
        L (float): Inductance in Henrys
        C (float): Capacitance in Farads
        Z0 (float): System impedance in Ohms (typically 50)
        frequency (float): Frequency in Hertz

    Returns:
        float: S11 value in dB
    """
    omega = 2 * np.pi * frequency
    # Impedance of the RLC circuit
    Z_coil = R + 1j * omega * L - 1j / (omega * C)
    # Reflection coefficient Gamma
    gamma = (Z_coil - Z0) / (Z_coil + Z0)
    # S11 in dB
    s11_db = 20 * np.log10(np.abs(gamma))
    return s11_db

# --- Simulation Parameters ---
# Let's model a proton resonance frequency for a 1.5T MRI scanner
f0 = 63.87e6  # Resonant frequency in Hz (~64 MHz)
Z0 = 50.0      # VNA characteristic impedance in Ohms

# --- Scenario 1: Properly Resonant Coil (High-Q) ---
# This represents a coil that is correctly biased or doesn't have active decoupling.
Q_high = 200.0      # Typical Q-factor for an unloaded receive coil
L = 250e-9          # Assume an inductance of 250 nH
# Calculate required C and R for the high-Q case
omega0 = 2 * np.pi * f0
C_high = 1 / (omega0**2 * L)
R_high = omega0 * L / Q_high

# --- Scenario 2: Decoupled Coil (Low-Q) ---
# This simulates an actively decoupled coil where a PIN diode adds high loss.
# The physical L and C are the same, but the effective resistance is much higher.
# Let's say the decoupling circuit adds ~20 Ohms of resistance.
R_low = R_high + 20.0
Q_low = omega0 * L / R_low

# --- Calculate and Print Results ---
s11_resonant_db = calculate_s11_db(R_high, L, C_high, Z0, f0)
s11_decoupled_db = calculate_s11_db(R_low, L, C_high, Z0, f0) # Note: use same C

print("--- Simulating VNA Measurement of an MRI Coil ---")
print(f"Target Resonance Frequency: {f0/1e6:.2f} MHz\n")

print("--- Scenario 1: Properly Resonant (High-Q) Coil ---")
print(f"Coil Resistance (R): {R_high:.4f} Ohms")
print(f"Coil Inductance (L): {L*1e9:.2f} nH")
print(f"Coil Capacitance (C): {C_high*1e12:.2f} pF")
print(f"Resulting Quality Factor (Q): {Q_high:.1f}")
print("Final Equation for S11 at resonance:")
print(f"S11({f0/1e6:.2f} MHz) = {s11_resonant_db:.2f} dB")
print("--> This shows a strong, deep resonance peak on the VNA.\n")


print("--- Scenario 2: Actively Decoupled (Low-Q) Coil ---")
print(f"Coil Resistance (R): {R_low:.4f} Ohms (Original R + decoupling resistance)")
print(f"Coil Inductance (L): {L*1e9:.2f} nH")
print(f"Coil Capacitance (C): {C_high*1e12:.2f} pF")
print(f"Resulting Quality Factor (Q): {Q_low:.1f}")
print("Final Equation for S11 at resonance:")
print(f"S11({f0/1e6:.2f} MHz) = {s11_decoupled_db:.2f} dB")
print("--> The resonance is almost completely gone, showing a nearly flat line on the VNA.")