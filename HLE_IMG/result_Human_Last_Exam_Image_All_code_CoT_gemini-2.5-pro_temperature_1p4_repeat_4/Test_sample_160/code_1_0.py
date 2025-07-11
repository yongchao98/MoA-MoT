import numpy as np

# 1. Define constants and parameters from the problem description and diagram
V_RF = 1.0          # V, peak voltage of the fundamental
f0 = 915e6          # Hz, fundamental frequency (from w = 2*pi*915MHz)
R0 = 50.0           # Ohms, base parasitic resistance
C_parasitic = 2e-15 # Farads (2 fF)
Cc = 1e-12          # Farads (1 pF)
R_L = 8000.0        # Ohms (8 kOhm)

# 2. State and apply simplifying assumptions
# Assume the effective AC input resistance of the rectifier is half the DC load resistance.
R_L_eff = R_L / 2.0
# Assume the total DC voltage is the sum of DC contributions from each harmonic,
# with a rectification gain of 1 (V_dc = V_peak_ac).

# 3. Calculate the source voltage for each harmonic
# The voltage drops by 10% for each higher harmonic relative to the previous one.
harmonics_source_voltages = {
    1: V_RF,
    3: V_RF * 0.9,
    5: V_RF * 0.9 * 0.9,
    7: V_RF * 0.9 * 0.9 * 0.9,
}

# 4. Calculate the attenuated voltage reaching the rectifier for each harmonic
rectified_voltages = {}
total_dc_voltage = 0

print("Calculating the final DC voltage by summing the contributions from each harmonic.")
print("-" * 75)

for n, v_source in harmonics_source_voltages.items():
    # Calculate frequency and angular frequency for the current harmonic
    f_n = n * f0
    w_n = 2 * np.pi * f_n

    # Calculate frequency-dependent parasitic resistance
    R_p_n = R0 * (f_n / f0)**2

    # Calculate impedance of the two series coupling capacitors (Cc)
    Z_Cc_n = 1 / (1j * w_n * Cc)
    Z_series_caps = 2 * Z_Cc_n

    # Total series impedance is the sum of parasitic resistance and coupling capacitor impedance
    Z_series_n = R_p_n + Z_series_caps

    # Calculate impedance of the shunt parasitic capacitance (C_parasitic)
    Z_Cp_n = 1 / (1j * w_n * C_parasitic)

    # Calculate the effective input impedance of the rectifier stage
    # This is the parallel combination of the effective load and the parasitic capacitance
    Z_rect_in_n = 1 / (1 / R_L_eff + 1 / Z_Cp_n)

    # Use the voltage divider formula to find the peak AC voltage at the rectifier input
    v_rect_n = v_source * np.abs(Z_rect_in_n / (Z_series_n + Z_rect_in_n))
    
    rectified_voltages[n] = v_rect_n
    total_dc_voltage += v_rect_n

# 5. Print the final results in the required equation format
v1 = rectified_voltages[1]
v3 = rectified_voltages[3]
v5 = rectified_voltages[5]
v7 = rectified_voltages[7]

print("\nFinal Calculation Summary:")
print("The total voltage across C_L is approximated by the sum of the rectified peak voltages:")
print(f"V_total = V_rect_1 + V_rect_3 + V_rect_5 + V_rect_7")
print(f"V_total = {v1:.4f}V + {v3:.4f}V + {v5:.4f}V + {v7:.4f}V")
print(f"V_total = {total_dc_voltage:.4f} V")

<<<2.7746>>>