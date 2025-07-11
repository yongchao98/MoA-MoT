import numpy as np
import math

# 1. Define all the given circuit parameters and constants.
V_RF = 1.0  # V, peak voltage of the fundamental
f0 = 915e6  # Hz, fundamental frequency
omega0 = 2 * np.pi * f0
Cc = 1e-12  # F, from table
CL = 5e-12  # F, from table
RL = 8e3   # Ohms, from table
R0_parasitic = 50.0  # Ohms, base parasitic resistance from text
f_ref_parasitic = 915e6 # Hz, reference frequency for parasitic resistance from text
C_parasitic = 2e-15 # F, from text

# 2. Define the harmonics to be considered.
harmonics = [1, 3, 5, 7]

# 3. Create a dictionary of input voltage amplitudes for each harmonic.
# The voltage drops by 10% for each higher harmonic relative to the previous harmonic.
V_in_amplitudes = {}
V_in_amplitudes[1] = V_RF
V_in_amplitudes[3] = V_in_amplitudes[1] * (1 - 0.1)
V_in_amplitudes[5] = V_in_amplitudes[3] * (1 - 0.1)
V_in_amplitudes[7] = V_in_amplitudes[5] * (1 - 0.1)

# 4. Loop through each harmonic to calculate its voltage contribution at the load.
load_voltage_peaks = []
# The effective load capacitance is the sum of the main load capacitor and the parasitic capacitor.
C_L_eff = CL + C_parasitic

print("Calculating the peak voltage at the load for each harmonic:")
print("-" * 50)

for n in harmonics:
    # a. Calculate the angular frequency for the current harmonic.
    omega_n = n * omega0

    # b. Calculate the frequency-dependent parasitic resistance.
    R_parasitic_n = R0_parasitic * (n * f0 / f_ref_parasitic)**2

    # c. Calculate the impedance of the coupling capacitor.
    Z_Cc_n = 1 / (1j * omega_n * Cc)

    # d. Calculate the total series impedance.
    Z_series_n = R_parasitic_n + Z_Cc_n

    # e. Calculate the impedance of the parallel RC load.
    # The admittance Y is the sum of the admittances of the parallel components. Z = 1/Y.
    Y_load_n = 1/RL + 1j * omega_n * C_L_eff
    Z_load_n = 1 / Y_load_n

    # f. Use the voltage divider formula to find the complex voltage across the load.
    V_in_n = V_in_amplitudes[n]
    V_load_n_complex = V_in_n * (Z_load_n / (Z_series_n + Z_load_n))

    # g. The peak voltage is the magnitude (absolute value) of the complex voltage.
    V_load_n_peak = np.abs(V_load_n_complex)
    load_voltage_peaks.append(V_load_n_peak)
    
    print(f"For harmonic n={n}:")
    print(f"  V_in = {V_in_n:.3f} V")
    print(f"  Z_series = {Z_series_n.real:.2f} + {Z_series_n.imag:.2f}j Ohms")
    print(f"  Z_load = {Z_load_n.real:.2f} + {Z_load_n.imag:.2f}j Ohms")
    print(f"  V_load_peak = {V_load_n_peak:.4f} V")
    print("-" * 50)

# 5. Sum the peak voltages of all harmonics to get the total voltage.
# This represents the total voltage contribution from all harmonic components.
total_voltage = sum(load_voltage_peaks)

# 6. Print the final result in the requested format.
print("The total voltage across C_L is approximated by the sum of the peak voltages:")
equation = "V_total = {:.4f} V + {:.4f} V + {:.4f} V + {:.4f} V = {:.4f} V".format(
    load_voltage_peaks[0],
    load_voltage_peaks[1],
    load_voltage_peaks[2],
    load_voltage_peaks[3],
    total_voltage
)
print(equation)

<<<0.1909>>>