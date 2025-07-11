import math

def calculate_resonant_frequency(L, C):
    """Calculates the resonant frequency of an LC circuit."""
    return 1 / (2 * math.pi * math.sqrt(L * C))

# --- Parameters for a typical 1.5T MRI coil ---
# Operating frequency for protons at 1.5 Tesla is ~63.8 MHz
L_coil = 150e-9  # Inductance of the coil in Henries (150 nH)
C_tune = 41.5e-12 # Tuning capacitance in Farads (41.5 pF)

# --- Parameters for the detuned state ---
# When the decoupling circuit is inactive (probed by a VNA),
# it often adds capacitance to the circuit, shifting the resonance down.
C_detune = 500e-12 # Detuning capacitance in Farads (500 pF)
C_total_detuned = C_tune + C_detune

# --- Calculations ---
# 1. Calculate the intended operating frequency (Tuned State)
freq_tuned = calculate_resonant_frequency(L_coil, C_tune)

# 2. Calculate the frequency in the detuned state (the state you are likely measuring)
freq_detuned = calculate_resonant_frequency(L_coil, C_total_detuned)

# --- Output the results ---
print("This script demonstrates why a functional MRI coil might not show resonance when probed.")
print("It calculates the resonant frequency in two states: 'Tuned' (active in scanner) and 'Detuned' (inactive, default state).\n")

print("--- Coil Parameters ---")
print(f"Coil Inductance (L): {L_coil * 1e9:.1f} nH")
print(f"Tuning Capacitance (C_tune): {C_tune * 1e12:.1f} pF")
print(f"Detuning Capacitance (C_detune): {C_detune * 1e12:.1f} pF\n")

print("--- Tuned State (Active in Scanner) ---")
print("This is the state when the coil receives a control signal from the MRI.")
print(f"f_res = 1 / (2 * pi * sqrt(L * C_tune))")
print(f"f_res = 1 / (2 * pi * sqrt({L_coil:.3e} H * {C_tune:.3e} F))")
print(f"Expected Resonant Frequency: {freq_tuned / 1e6:.2f} MHz\n")

print("--- Detuned State (Inactive, Probed by VNA) ---")
print("This is the default state without the MRI control signal.")
print(f"C_total = C_tune + C_detune = {C_total_detuned * 1e12:.1f} pF")
print(f"f_res = 1 / (2 * pi * sqrt(L * C_total))")
print(f"f_res = 1 / (2 * pi * sqrt({L_coil:.3e} H * {C_total_detuned:.3e} F))")
print(f"Actual Measured Resonant Frequency: {freq_detuned / 1e6:.2f} MHz\n")

print("Conclusion:")
print("Your VNA is likely measuring the coil in its 'Detuned State',")
print("where the resonance is shifted to a much lower frequency, far from the expected operating frequency.")
