import cmath
import math

# --- 1. Define Constants ---
V_RF = 1.0  # V, peak amplitude of the fundamental
f1 = 915e6  # Hz, fundamental frequency
C_C = 1e-12  # F, coupling capacitance
R_L = 8000  # Ohms, load resistance
# C_L = 5e-12 F (This component's primary role is DC filtering and is implicitly included in the R_L/2 load approximation)
R0 = 50.0  # Ohms, parasitic resistance base value
f0 = 915e6  # Hz, parasitic resistance reference frequency
C_parasitic = 2e-15  # F, parasitic capacitance

# --- 2. Calculate Losses at the Fundamental Frequency (f1) ---
w1 = 2 * math.pi * f1

# Parasitic resistance at f1
R_p_f1 = R0 * (f1 / f0)**2

# Impedance of the two series coupling capacitors
Z_Cc_f1 = 1 / (1j * w1 * C_C)
Z_Cc_total_f1 = 2 * Z_Cc_f1

# The effective input resistance of the rectifier is approximated as R_L/2
R_eq = R_L / 2

# Impedance of the parasitic capacitance at f1
Z_Cpara_f1 = 1 / (1j * w1 * C_parasitic)

# The total series impedance in the signal path
Z_series = R_p_f1 + Z_Cc_total_f1

# The input impedance of the rectifier stage (load of the divider)
Z_in_rect = 1 / (1/R_eq + 1/Z_Cpara_f1)

# The total impedance seen by the source
Z_total = Z_series + Z_in_rect

# Calculate the voltage divider gain
H_f1 = Z_in_rect / Z_total
gain = abs(H_f1)

# --- 3. Calculate Final DC Voltage ---

# Calculate the DC voltage from the fundamental component, considering circuit losses
V_DC_fundamental = V_RF * gain * 2  # Factor of 2 for the voltage doubler

# Define harmonic loss factors
loss_factor_3rd = 0.9  # 1 - 10%
loss_factor_5th = 0.9  # 1 - 10%
loss_factor_7th = 0.9  # 1 - 10%

# Apply harmonic distortion losses
V_DC_final = V_DC_fundamental * loss_factor_3rd * loss_factor_5th * loss_factor_7th

# --- 4. Print the Calculation Steps and Final Result ---

print("Calculation of the total voltage across the load capacitor CL:")
print("-" * 60)

print("Step 1: Calculate the gain from circuit and parasitic losses at the fundamental frequency (915 MHz).")
print(f"The gain |H(f1)| is calculated as |Z_in_rect / (Z_series + Z_in_rect)|.")
print(f"Gain |H(f1)| = {gain:.4f}\n")

print("Step 2: Calculate the ideal DC voltage from the fundamental frequency, including circuit losses.")
print("V_DC_fundamental = V_RF * gain * 2 (for voltage doubler)")
print(f"V_DC_fundamental = {V_RF:.2f} V * {gain:.4f} * 2 = {V_DC_fundamental:.4f} V\n")

print("Step 3: Apply successive voltage drops due to harmonic distortions.")
print("V_DC_final = V_DC_fundamental * (1 - 0.1) * (1 - 0.1) * (1 - 0.1)")
final_equation_str = (
    f"V_DC_final = {V_DC_fundamental:.4f} V * "
    f"{loss_factor_3rd} * {loss_factor_5th} * {loss_factor_7th}"
)
print(final_equation_str)
print(f"V_DC_final = {V_DC_final:.4f} V\n")

print("The total voltage across C_L is approximately:")
print(f"{V_DC_final:.3f} V")
print("-" * 60)
print(f'<<<1.429>>>')