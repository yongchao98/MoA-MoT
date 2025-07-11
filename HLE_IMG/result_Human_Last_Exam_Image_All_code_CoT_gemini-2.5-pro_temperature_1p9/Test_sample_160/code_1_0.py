import math

# Step 1: Define all the given constants from the problem description and diagram.
V_rf_peak_fundamental = 1.0  # V
f_fundamental = 915e6  # Hz, from omega = 2*pi*915MHz
R_L = 8e3  # Ohms
C_L = 5e-12  # Farads
R0_parasitic = 50.0  # Ohms
f0_parasitic = 915e6  # Hz
C_parasitic = 2e-15  # Farads

# Assumptions based on the problem and standard RF analysis
# The rectifier's input impedance is assumed to be 50 Ohms, a common design target.
R_rectifier_input = 50.0 # Ohms
# The cross-coupled rectifier circuit acts as a voltage doubler.
voltage_gain = 2.0

print("Step-by-step Calculation of the Total Voltage Across C_L\n")

# Step 2: Calculate the peak voltage for the fundamental and each harmonic.
print("--- Step 1: Harmonic Input Voltages ---")
harmonics = [1, 3, 5, 7]
V_peaks = {}
V_peaks[1] = V_rf_peak_fundamental
print(f"Fundamental (1st) harmonic peak voltage: {V_peaks[1]:.3f} V")
for i in range(1, len(harmonics)):
    prev_harmonic = harmonics[i-1]
    current_harmonic = harmonics[i]
    # Voltage drops by 10% relative to the previous harmonic
    V_peaks[current_harmonic] = V_peaks[prev_harmonic] * 0.9
    print(f"{current_harmonic}rd/th harmonic peak voltage: {V_peaks[current_harmonic]:.3f} V")

# Step 3: Calculate the frequency-dependent parasitic resistance for each harmonic.
print("\n--- Step 2: Frequency-Dependent Parasitic Resistances ---")
R_parasitics = {}
frequencies = {}
for n in harmonics:
    f_n = n * f_fundamental
    frequencies[n] = f_n
    # R_parasitic(f) = R0 * (f/f0)^2
    R_parasitics[n] = R0_parasitic * (f_n / f0_parasitic)**2
    print(f"Parasitic resistance at f_{n} ({n*f_fundamental/1e6:.0f} MHz): {R_parasitics[n]:.1f} Ohms")

# Step 4: Calculate the attenuated voltage at the rectifier input for each harmonic.
print("\n--- Step 3: Attenuated Voltages at Rectifier Input (Voltage Divider) ---")
V_attenuated_peaks = {}
for n in harmonics:
    # V_attenuated = V_peak * (R_rect_in / (R_rect_in + R_parasitic))
    V_attenuated_peaks[n] = V_peaks[n] * (R_rectifier_input / (R_rectifier_input + R_parasitics[n]))
    print(f"Attenuated peak voltage for {n}st/rd/th harmonic: {V_attenuated_peaks[n]:.4f} V")

# Step 5: Sum the attenuated peak voltages and apply the rectifier's voltage gain.
print("\n--- Step 4 & 5: Final Output Voltage Calculation ---")
total_attenuated_peak_voltage = sum(V_attenuated_peaks.values())
print(f"Sum of attenuated peak voltages = ", " + ".join([f"{v:.4f}" for v in V_attenuated_peaks.values()]), f" = {total_attenuated_peak_voltage:.4f} V")

# The rectifier acts as a voltage doubler.
# The load (R_L, C_L, C_parasitic) has a long time constant, so V_DC is close to the peak.
final_voltage_across_CL = total_attenuated_peak_voltage * voltage_gain

print(f"\nFinal calculated voltage across C_L (V_DC) = Sum of attenuated voltages * Voltage Gain")
print(f"V_DC = {total_attenuated_peak_voltage:.4f} V * {voltage_gain}")
print(f"Total Voltage Across C_L = {final_voltage_across_CL:.4f} V")

# Final answer in the required format
# <<<1.2715>>>