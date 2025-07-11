import math

# --- Step 1: Define constants and parameters from the problem ---
V_RF_fundamental = 1.0  # V, amplitude of the fundamental
f0 = 915e6  # Hz, fundamental frequency
R0_parasitic = 50.0  # Ohms, parasitic resistance base value
f0_parasitic_ref = 915e6  # Hz, parasitic resistance reference frequency
C_parasitic = 2e-15  # F, parasitic capacitance (2 fF)
R_L = 8000.0  # Ohms, load resistance
# Assumption: The rectifier input impedance is matched to a standard 50 Ohm system
Z_rectifier_assumed = 50.0 # Ohms

print("--- System Parameters ---")
print(f"Fundamental Voltage (V1): {V_RF_fundamental} V")
print(f"Fundamental Frequency (f0): {f0/1e6} MHz")
print(f"Parasitic Resistance (R0): {R0_parasitic} Ohms at {f0_parasitic_ref/1e6} MHz")
print(f"Parasitic Capacitance (Cp): {C_parasitic*1e15} fF")
print(f"Assumed Rectifier Impedance (Z_rect): {Z_rectifier_assumed} Ohms\n")


# --- Step 2: Calculate Harmonic Content Efficiency (eta_h) ---
print("--- Step 2: Harmonic Content Efficiency Calculation ---")
# Calculate harmonic voltage amplitudes
V1 = V_RF_fundamental
V3 = V1 * 0.9
V5 = V3 * 0.9
V7 = V5 * 0.9
print(f"Harmonic Voltages: V3 = {V3:.3f} V, V5 = {V5:.3f} V, V7 = {V7:.3f} V")

# Calculate power ratios (proportional to V^2)
P1_rel = V1**2
P3_rel = V3**2
P5_rel = V5**2
P7_rel = V7**2

P_total_rel = P1_rel + P3_rel + P5_rel + P7_rel
eta_h = P1_rel / P_total_rel

print(f"Relative Power of Fundamental (V1^2): {P1_rel:.4f}")
print(f"Total Relative Power (V1^2 + V3^2 + V5^2 + V7^2): {P1_rel:.4f} + {P3_rel:.4f} + {P5_rel:.4f} + {P7_rel:.4f} = {P_total_rel:.4f}")
print(f"Harmonic Efficiency Factor (eta_h) = P1_rel / P_total_rel = {eta_h:.4f}\n")


# --- Step 3: Calculate Parasitic Transfer Efficiency (eta_p) for the fundamental ---
print("--- Step 3: Parasitic Transfer Efficiency Calculation (at fundamental frequency) ---")
f1 = f0
# Calculate parasitic resistance at f1
R_p1 = R0_parasitic * (f1 / f0_parasitic_ref)**2
print(f"Parasitic Resistance at {f1/1e6} MHz (Rp1): {R_p1:.4f} Ohms")

# Calculate impedance of the parallel load (Z_rectifier || C_parasitic)
omega1 = 2 * math.pi * f1
# Admittance of C_parasitic is j*omega1*C_parasitic
# Admittance of Z_rectifier is 1/Z_rectifier_assumed
# Total load admittance Y_L = 1/Z_rectifier + j*omega1*C_parasitic
# Total load impedance Z_L = 1 / Y_L
# Power transfer efficiency depends on the real parts of impedances
# PTE = Re(Z_L) / (R_p1 + Re(Z_L))
# Re(Z_L) = Re(1 / (1/R_rect + j*w*C)) = R_rect / (1 + (w*C*R_rect)^2)
Re_Z_L1 = Z_rectifier_assumed / (1 + (omega1 * C_parasitic * Z_rectifier_assumed)**2)
print(f"Effective Load Resistance at {f1/1e6} MHz (Re(Z_L1)): {Re_Z_L1:.4f} Ohms")

# Calculate power transfer efficiency
eta_p_f1 = Re_Z_L1 / (R_p1 + Re_Z_L1)
print(f"Parasitic Transfer Efficiency (eta_p_f1) = Re(Z_L1) / (Rp1 + Re(Z_L1)) = {eta_p_f1:.4f}\n")


# --- Step 4: Calculate Overall System Efficiency ---
print("--- Step 4: Overall System Efficiency Calculation ---")
# Assuming intrinsic rectifier efficiency is 1 (ideal conversion of received fundamental power)
eta_intrinsic = 1.0
eta_system = eta_intrinsic * eta_p_f1 * eta_h

print(f"Overall Efficiency (eta_system) = eta_intrinsic * eta_p_f1 * eta_h")
print(f"eta_system = {eta_intrinsic:.2f} * {eta_p_f1:.4f} * {eta_h:.4f}")
print(f"eta_system = {eta_system:.4f}")
print(f"Overall System Efficiency: {eta_system * 100:.2f}%")

# Final Answer Formatting
final_answer = eta_system * 100
# print(f'<<<{final_answer:.2f}>>>')