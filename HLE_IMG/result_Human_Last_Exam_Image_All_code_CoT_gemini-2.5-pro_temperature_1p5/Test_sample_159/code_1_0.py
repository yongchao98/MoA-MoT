import math

# --- 1. Define Constants and Input Parameters ---
V1_amp = 1.0  # V, amplitude of the fundamental frequency
voltage_drop_factor = 0.10 # 10% drop for each higher harmonic
R0 = 50.0 # Ohms, characteristic parasitic resistance
f0 = 915e6 # Hz, reference frequency for parasitic resistance
# Assumption: Rectifier is matched for max power transfer at the fundamental frequency.
R_rectifier_assumed = 50.0 # Ohms

# --- 2. Calculate Harmonic Voltages ---
print("Step 1: Calculate harmonic voltage amplitudes.")
V3_amp = V1_amp * (1 - voltage_drop_factor)
V5_amp = V3_amp * (1 - voltage_drop_factor)
V7_amp = V5_amp * (1 - voltage_drop_factor)

print(f"V1 amplitude = {V1_amp:.3f} V")
print(f"V3 amplitude = {V1_amp:.3f} * (1 - {voltage_drop_factor}) = {V3_amp:.3f} V")
print(f"V5 amplitude = {V3_amp:.3f} * (1 - {voltage_drop_factor}) = {V5_amp:.3f} V")
print(f"V7 amplitude = {V5_amp:.3f} * (1 - {voltage_drop_factor}) = {V7_amp:.3f} V\n")

# --- 3. Calculate Power Proportions and Harmonic Efficiency (eta_h) ---
print("Step 2: Calculate harmonic efficiency (eta_h).")
print("Assuming constant input impedance for all harmonics, Power is proportional to Voltage^2.")
# Power is proportional to V^2
P1_prop = V1_amp**2
P3_prop = V3_amp**2
P5_prop = V5_amp**2
P7_prop = V7_amp**2

P_total_prop = P1_prop + P3_prop + P5_prop + P7_prop

print(f"Relative power of fundamental (P1): {V1_amp:.3f}^2 = {P1_prop:.4f}")
print(f"Relative power of 3rd harmonic (P3): {V3_amp:.3f}^2 = {P3_prop:.4f}")
print(f"Relative power of 5th harmonic (P5): {V5_amp:.3f}^2 = {P5_prop:.4f}")
print(f"Relative power of 7th harmonic (P7): {V7_amp:.3f}^2 = {P7_prop:.4f}")
print(f"Total relative power (P_total): {P1_prop:.4f} + {P3_prop:.4f} + {P5_prop:.4f} + {P7_prop:.4f} = {P_total_prop:.4f}\n")

# eta_h is the ratio of useful power (fundamental) to total power
eta_h = P1_prop / P_total_prop
print("Harmonic efficiency assumes only fundamental power is converted.")
print(f"eta_h = P1 / P_total = {P1_prop:.4f} / {P_total_prop:.4f} = {eta_h:.4f}\n")


# --- 4. Calculate Parasitic Resistance and Efficiency (eta_p) ---
print("Step 3: Calculate parasitic loss efficiency (eta_p) at the fundamental frequency.")
f1 = 915e6 # Hz, fundamental frequency
R_parasitic_f1 = R0 * (f1 / f0)**2
print(f"Parasitic resistance at f1={f1/1e6} MHz: R_p(f1) = {R0} * ({f1/1e6}/{f0/1e6})^2 = {R_parasitic_f1:.2f} Ohms")

print(f"Assuming rectifier input resistance R_rectifier = {R_rectifier_assumed:.2f} Ohms (for matching).")
# eta_p is the efficiency of the power divider formed by R_parasitic and R_rectifier
eta_p = R_rectifier_assumed / (R_parasitic_f1 + R_rectifier_assumed)
print("Parasitic efficiency is the fraction of power delivered to the rectifier.")
print(f"eta_p = R_rectifier / (R_p(f1) + R_rectifier) = {R_rectifier_assumed:.2f} / ({R_parasitic_f1:.2f} + {R_rectifier_assumed:.2f}) = {eta_p:.4f}\n")

# --- 5. Calculate Overall System Efficiency ---
print("Step 4: Calculate the overall system efficiency.")
eta_total = eta_h * eta_p
print("Total Efficiency = Harmonic Efficiency * Parasitic Loss Efficiency")
print(f"eta_total = eta_h * eta_p = {eta_h:.4f} * {eta_p:.4f} = {eta_total:.4f}")

# Express as a percentage
eta_total_percent = eta_total * 100
print(f"\nThe overall system efficiency is {eta_total_percent:.2f}%.")
print(f"\nFinal calculation: Overall Efficiency = ({P1_prop:.4f} / {P_total_prop:.4f}) * ({R_rectifier_assumed:.2f} / ({R_parasitic_f1:.2f} + {R_rectifier_assumed:.2f})) = {eta_total:.4f}")
print(f"<<<{eta_total:.4f}>>>")