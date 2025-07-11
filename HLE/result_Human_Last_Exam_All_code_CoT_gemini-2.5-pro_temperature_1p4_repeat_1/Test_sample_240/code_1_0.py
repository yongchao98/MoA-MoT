import math

def calculate_resonant_frequency(inductance_H, capacitance_F):
    """Calculates the resonant frequency of an LC circuit in MHz."""
    if inductance_H <= 0 or capacitance_F <= 0:
        return 0
    # Formula for resonance: f = 1 / (2 * pi * sqrt(L * C))
    frequency_hz = 1 / (2 * math.pi * math.sqrt(inductance_H * capacitance_F))
    frequency_mhz = frequency_hz / 1e6
    return frequency_mhz

# --- State 1: Tuned Coil (in "Receive" mode) ---
# Let's assume our coil is for a 1.5T MRI system (~63.87 MHz)
# We can pick a realistic inductance, e.g., 150 nanohenries (nH)
L_tuned_H = 150e-9  # 150 nH in Henries

# We can calculate the required capacitance to achieve the target frequency
# C = 1 / ( (2*pi*f)^2 * L )
target_freq_hz = 63.87e6
C_tuned_F = 1 / ((2 * math.pi * target_freq_hz)**2 * L_tuned_H) # Should be around 41.3 pF

# Let's use a clean number for the demonstration
C_tuned_F = 41.3e-12 # 41.3 picoFarads in Farads

# Calculate the actual frequency
resonant_freq_mhz = calculate_resonant_frequency(L_tuned_H, C_tuned_F)

print("--- Step 1: Simulating a normal, tuned MRI coil ---")
print("The coil requires a specific DC bias to be in this 'receive' mode.")
print(f"Inductance (L): {L_tuned_H * 1e9:.1f} nH")
print(f"Capacitance (C): {C_tuned_F * 1e12:.1f} pF")
print("\nCalculating the resonant frequency:")
print(f"f = 1 / (2 * pi * sqrt(L * C))")
print(f"f = 1 / (2 * {math.pi:.4f} * sqrt({L_tuned_H:.3e} H * {C_tuned_F:.3e} F))")
print(f"Result: The coil is resonant at {resonant_freq_mhz:.2f} MHz. This peak would be visible on a VNA if the coil were active.\n")


# --- State 2: Decoupled Coil (Default state without DC bias) ---
# A PIN diode circuit effectively changes the capacitance. Let's simulate this
# by putting a tiny 1 pF capacitor in series with our tuning capacitor.
# The new total capacitance C_new = (C1 * C2) / (C1 + C2)
C_pin_diode_effect_F = 1.0e-12 # 1.0 pF
C_decoupled_F = (C_tuned_F * C_pin_diode_effect_F) / (C_tuned_F + C_pin_diode_effect_F)

# Calculate the new, decoupled frequency
decoupled_freq_mhz = calculate_resonant_frequency(L_tuned_H, C_decoupled_F)

print("--- Step 2: Simulating the decoupled coil (no DC bias) ---")
print("Without the DC bias, the PIN diode circuit alters the capacitance.")
print(f"Original Capacitance: {C_tuned_F * 1e12:.1f} pF")
print(f"New Effective Capacitance (C_decoupled): {C_decoupled_F * 1e12:.3f} pF")
print("\nCalculating the new 'decoupled' resonant frequency:")
print(f"f_decoupled = 1 / (2 * pi * sqrt(L * C_decoupled))")
print(f"f_decoupled = 1 / (2 * {math.pi:.4f} * sqrt({L_tuned_H:.3e} H * {C_decoupled_F:.3e} F))")
print(f"Result: The resonance shifts to {decoupled_freq_mhz:.2f} MHz.\n")

print("--- Conclusion ---")
print(f"The resonant frequency has been shifted from {resonant_freq_mhz:.2f} MHz to {decoupled_freq_mhz:.2f} MHz.")
print("This new frequency is far outside the expected range and would not be seen on the VNA, making the coil appear non-resonant.")
