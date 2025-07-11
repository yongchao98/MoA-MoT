import cmath

# --- Preamplifier Characteristics ---
# These are typical S-parameters for a Low Noise Amplifier (LNA)
# S11: Input reflection coefficient (how much is reflected from the input)
# S21: Forward gain (how much the signal is amplified going forward)
# S12: Reverse isolation (how much the signal is attenuated going backward)
# S22: Output reflection coefficient (how much is reflected from the output)

# Let's define them in decibels (dB), which is common
S11_dB = -15  # Good input match
S21_dB = 20   # 20 dB of gain
S12_dB = -45  # Very high reverse isolation, this is the key parameter!
S22_dB = -12  # Good output match

# --- VNA Test Signal ---
# Let's assume the VNA sends a signal with a power of 0 dBm.
# In linear magnitude, 0 dBm is 1. A VNA measures the ratio of signals.
vna_signal_magnitude = 1.0

# --- Function to convert dB to linear magnitude ---
def db_to_linear(db_value):
    """Converts a value in decibels (dB) to a linear magnitude ratio."""
    return 10**(db_value / 20.0)

# --- Calculation ---
# Convert the dB S-parameters to linear magnitude ratios
s21_linear = db_to_linear(S21_dB)
s12_linear = db_to_linear(S12_dB)

# Scenario 1: A real MR signal from the coil goes FORWARD through the preamp
# The input is the tiny MR signal, let's say its magnitude is 'x'
# The output signal would be x * s21_linear
# Let's calculate the gain factor
forward_gain_factor = s21_linear

# Scenario 2: The VNA signal tries to go BACKWARD through the preamp
# The input is the VNA signal from the connector, let's say its magnitude is 1.0
# The signal reaching the resonant element would be 1.0 * s12_linear
# Let's calculate the attenuation factor
backward_attenuation_factor = s12_linear

# --- Output the results ---
print("--- Preamplifier Signal Path Analysis ---")
print(f"A signal going FORWARD (Coil -> Connector) is multiplied by a factor of: {forward_gain_factor:.2f}")
print(f"This is an amplification of {S21_dB} dB.")
print("\n")
print(f"A signal going BACKWARD (VNA -> Coil) is multiplied by a factor of: {backward_attenuation_factor:.8f}")
print(f"This is an attenuation of {abs(S12_dB)} dB.")
print("\n--- Conclusion ---")
print("The VNA's signal is so heavily attenuated by the preamplifier's reverse isolation")
print("that it's too weak to interact with the resonant element of the coil.")
print("Therefore, the VNA cannot 'see' the coil's resonance.")
