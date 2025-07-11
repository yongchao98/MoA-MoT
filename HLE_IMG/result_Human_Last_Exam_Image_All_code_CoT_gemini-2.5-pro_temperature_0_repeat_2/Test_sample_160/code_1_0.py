import numpy as np

# --- 1. Define Constants from the problem ---
V_RF_fund = 1.0  # V, amplitude of the fundamental input
f_fund = 915e6   # Hz, fundamental frequency
Cc = 1e-12       # F, coupling capacitance
CL = 5e-12       # F, load capacitance
RL = 8000        # Ohms, load resistance
R0_parasitic = 50 # Ohms, base parasitic resistance
f0_parasitic = 915e6 # Hz, reference frequency for parasitic resistance
# C_parasitic = 2e-15 F is negligible

harmonics = [1, 3, 5, 7]
V_rect_in_n_list = []
phi_n_list = []

print("--- Step-by-Step Calculation ---")
print("Calculating attenuated voltage for each harmonic...\n")

# --- 2. Calculate attenuated voltage for each harmonic ---
for n in harmonics:
    f_n = n * f_fund
    w_n = 2 * np.pi * f_n

    # Calculate harmonic input voltage amplitude
    V_rf_n = V_RF_fund * (0.9)**((n - 1) / 2)

    # Calculate frequency-dependent parasitic resistance
    R_p_n = R0_parasitic * (f_n / f0_parasitic)**2

    # Calculate magnitude of the impedance of the two series coupling capacitors
    Z_cc_n_mag = 2 / (w_n * Cc)

    # Calculate attenuated voltage amplitude at rectifier input (voltage divider)
    V_rect_in_n = V_rf_n * Z_cc_n_mag / np.sqrt(R_p_n**2 + Z_cc_n_mag**2)
    V_rect_in_n_list.append(V_rect_in_n)

    # Calculate phase shift in radians
    phi_n = np.arctan(Z_cc_n_mag / R_p_n) - np.pi / 2
    phi_n_list.append(phi_n)

    print(f"For Harmonic n={n}:")
    print(f"  - Input Voltage (V_rf_{n}): {V_rf_n:.3f} V")
    print(f"  - Parasitic Resistance (R_p_{n}): {R_p_n:.1f} Ohms")
    print(f"  - Rectifier Input Impedance (|Z_rect_{n}|): {Z_cc_n_mag:.1f} Ohms")
    print(f"  - Attenuated Voltage (V'_n): {V_rect_in_n:.4f} V")
    print(f"  - Phase Shift (phi_{n}): {np.degrees(phi_n):.1f} degrees\n")

# --- 3. Numerically find the peak of the composite waveform ---
T_fund = 1 / f_fund
# Use a high number of points for accuracy
time_points = np.linspace(0, T_fund, 5000)
v_rect_in_t = np.zeros_like(time_points)

for i, n in enumerate(harmonics):
    w_n = 2 * np.pi * n * f_fund
    v_rect_in_t += V_rect_in_n_list[i] * np.sin(w_n * time_points + phi_n_list[i])

V_peak_eff = np.max(v_rect_in_t)

print("--- Effective Peak Voltage Calculation ---")
print(f"The composite waveform at the rectifier input is the sum of the attenuated harmonics.")
print(f"By sampling this waveform, the effective peak voltage is found.")
print(f"Effective Peak Voltage (V_peak_eff): {V_peak_eff:.4f} V\n")


# --- 4. Calculate the final DC voltage ---
# Using the formula for a loaded voltage doubler
# V_DC = (2 * V_peak_eff) / (1 + 1 / (2 * f * C * R_L))
denominator_term = 1 / (2 * f_fund * Cc * RL)
V_DC = (2 * V_peak_eff) / (1 + denominator_term)

print("--- Final DC Voltage Calculation ---")
print("The DC voltage is calculated using the loaded voltage doubler formula:")
print(f"V_DC = (2 * V_peak_eff) / (1 + 1 / (2 * f_1 * Cc * R_L))")
print(f"V_DC = (2 * {V_peak_eff:.4f}) / (1 + 1 / (2 * {f_fund:.0f} * {Cc} * {RL}))")
print(f"V_DC = {2 * V_peak_eff:.4f} / (1 + {denominator_term:.4f})")
print(f"V_DC = {2 * V_peak_eff:.4f} / {1 + denominator_term:.4f}")
print(f"Total Voltage across C_L (V_DC) = {V_DC:.4f} V")

# Final answer tag
# print(f"\n<<<{V_DC:.4f}>>>")