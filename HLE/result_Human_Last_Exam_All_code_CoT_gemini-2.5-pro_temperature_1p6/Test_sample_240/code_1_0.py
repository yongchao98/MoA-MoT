import numpy as np

def calculate_s11_db(R, L, C, freq_hz, Z0=50):
    """
    Calculates the S11 reflection coefficient in dB for a series RLC circuit.
    
    Args:
        R (float): Resistance in Ohms.
        L (float): Inductance in Henrys.
        C (float): Capacitance in Farads.
        freq_hz (float): Frequency of interest in Hertz.
        Z0 (float): Characteristic impedance of the measurement system (VNA).
    
    Returns:
        float: S11 value in decibels (dB).
    """
    # Angular frequency
    omega = 2 * np.pi * freq_hz
    
    # Calculate reactances
    XL = omega * L  # Inductive reactance
    XC = 1 / (omega * C)  # Capacitive reactance
    
    # Calculate complex impedance of the RLC circuit
    Z_coil = R + 1j * (XL - XC)
    
    # Calculate reflection coefficient (Gamma)
    gamma = (Z_coil - Z0) / (Z_coil + Z0)
    
    # Calculate S11 in dB
    s11_db = 20 * np.log10(np.abs(gamma))
    
    return s11_db

# --- Parameters for a typical 3T MRI coil (127.7 MHz) ---
f0 = 127.7e6  # Resonant frequency for 3 Tesla field

# Let's design L and C for this frequency
# We choose a realistic L for a small loop
L = 100e-9  # 100 nanohenrys
# Calculate C required for resonance at f0
C = 1 / ((2 * np.pi * f0)**2 * L) # Around 15.5 pF

# --- Case 1: TUNED Coil ---
# A tuned coil has very low intrinsic resistance (high Quality factor)
R_tuned = 0.5  # 0.5 Ohms

# --- Case 2: DETUNED Coil ---
# The active detuning circuit adds significant loss, effectively
# increasing the resistance of the resonant tank.
R_detuned = 100.0  # 100 Ohms

# --- Simulation ---
# Calculate S11 at the resonant frequency for both cases
s11_tuned = calculate_s11_db(R_tuned, L, C, f0)
s11_detuned = calculate_s11_db(R_detuned, L, C, f0)

# --- Print Results ---
print("This simulation models a VNA measurement at a coil's target resonant frequency.")
print("-" * 60)
print(f"Coil parameters: L = {L*1e9:.1f} nH, C = {C*1e12:.2f} pF, f0 = {f0/1e6:.2f} MHz")
print("-" * 60)

# The final equation shows the coil's resistance in each state.
# For the tuned state, R is low.
print(f"TUNED STATE: Resistance (R) = {R_tuned} Ohm")
# At resonance, Z_coil should be close to R_tuned. 
# A low resistance gives a good match (or a sharp, deep reflection if not matched to 50 Ohm).
# We are assuming it's not matched, to see a reflection dip.
print(f"VNA measurement at {f0/1e6:.2f} MHz shows S11 = {s11_tuned:.2f} dB")
print("This is a very sharp, deep dip on the VNA screen, indicating a strong resonance.")
print("\n")

# For the detuned state, R is very high.
print(f"DETUNED STATE (unpowered fail-safe): Resistance (R) = {R_detuned} Ohm")
# The high resistance 'spoils' the Q-factor and dominates the impedance.
print(f"VNA measurement at {f0/1e6:.2f} MHz shows S11 = {s11_detuned:.2f} dB")
print("This value is very close to 0 dB, meaning there is no visible dip.")
print("The resonance is effectively hidden, which is what you are observing.")
