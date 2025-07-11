import math

# Step 1 & 2: Analyze Input Signal and Calculate Total Input Power
print("--- Step 1 & 2: Input Power Analysis ---")
# Given parameters from the problem description and the image table
V_RF_1 = 1.0  # V, amplitude of the fundamental input signal
R_s = 50.0 # Ohms, assumed source impedance, consistent with R0

# The voltage drops by 10% for each higher harmonic relative to the previous.
# This is interpreted as V_n = V_{n-2} * (1 - 0.10)
V = {}
V[1] = V_RF_1
V[3] = V[1] * 0.9
V[5] = V[3] * 0.9
V[7] = V[5] * 0.9

harmonics = [1, 3, 5, 7]
P_avail = {}
P_avail_total = 0.0

# Available power from a source is P_n = V_n^2 / (2*R_s)
for n in harmonics:
    P_avail[n] = V[n]**2 / (2 * R_s)
    P_avail_total += P_avail[n]
    print(f"Voltage of {n}-th harmonic V_{n}: {V[n]:.3f} V")
    print(f"Available power of {n}-th harmonic P_avail_{n}: {P_avail[n]*1000:.3f} mW")
print(f"Total available input power P_in_total: {P_avail_total*1000:.4f} mW\n")

# Step 3 & 4: Model Parasitics and Calculate Power Delivered to Rectifier
print("--- Step 3 & 4: Parasitic Loss Analysis ---")
f0 = 915e6  # Hz, reference frequency from the problem
R0 = 50.0  # Ohms, reference parasitic resistance from the problem
C_parasitic = 2e-15  # F, parasitic capacitance from the problem
Z_in = 50.0 # Ohms, assumed rectifier input impedance for matching at f1

# Calculate parasitic resistance at the fundamental frequency (f1 = f0)
f1 = f0
R_p_f1 = R0 * (f1 / f0)**2

# The reactance of the 2fF parasitic capacitor at 915MHz is extremely high (~87 kOhms).
# Assuming it's in shunt, its effect is negligible. We proceed by considering only the series parasitic resistance.
# The effective load seen by the source is R_load_eff = R_p_f1 + Z_in
R_load_eff = R_p_f1 + Z_in

# Calculate power transferred from the source to the combined load (R_p + Z_in)
P_del_to_load_eff = P_avail[1] * (4 * R_s * R_load_eff) / (R_s + R_load_eff)**2
# The power delivered to Z_in is a fraction of this, based on the voltage divider R_p vs Z_in
P_1_in_rect = P_del_to_load_eff * (Z_in / (R_p_f1 + Z_in))
print(f"Parasitic resistance R_p at f1={f1/1e6} MHz: {R_p_f1:.1f} Ohms")
print(f"Power at fundamental delivered to rectifier input P_1_in_rect: {P_1_in_rect*1000:.4f} mW\n")

# Step 5: Calculate DC Output Power
print("--- Step 5: Output Power Calculation ---")
R_L = 8000.0  # Ohms, load resistance from the table

# Power at rectifier input: P_1_in_rect = (V_in_se_peak)^2 / (2 * Z_in)
V_in_se_peak = math.sqrt(2 * Z_in * P_1_in_rect)
# The rectifier is a differential voltage doubler.
V_in_diff_peak = 2 * V_in_se_peak
# For an ideal rectifier, DC output voltage equals the peak differential input voltage.
V_DC = V_in_diff_peak
# DC output power is P_out = V_DC^2 / R_L
P_out = V_DC**2 / R_L
print(f"Peak differential voltage at rectifier input V_in_diff_peak: {V_in_diff_peak:.4f} V")
print(f"Ideal DC output voltage V_DC: {V_DC:.4f} V")
print(f"DC output power P_out: {P_out*1000:.4f} mW\n")

# Step 6: Calculate Overall Efficiency
print("--- Step 6: Overall System Efficiency ---")
efficiency = P_out / P_avail_total
print("Final Equation: efficiency = P_out / P_in_total")
print(f"efficiency = {P_out:.6f} W / {P_avail_total:.6f} W")
print(f"Final Efficiency = {efficiency:.6f}")
print(f"Overall system efficiency (%): {efficiency * 100:.4f} %")

# Final answer in the required format
final_answer = efficiency
# <<<0.007413>>>