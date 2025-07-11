import cmath
import math

# Step 1: Define constants from the problem description and image
V_RF_fundamental = 1.0  # V
f0 = 915e6  # Hz
R0 = 50.0  # Ohms
C_parasitic = 2e-15  # Farads (2 fF)
Cc = 1e-12  # Farads (1 pF)
CL = 5e-12  # Farads (5 pF)
RL = 8e3 # Ohms (Note: RL is high and its impedance will be ignored in AC analysis compared to CL)

# Step 2: Calculate source harmonic amplitudes
harmonics = [1, 3, 5, 7]
V_source = {}
V_current = V_RF_fundamental
print("Calculating source harmonic voltages:")
for n in harmonics:
    if n > 1:
        # Voltage drops by 10% for each higher harmonic
        V_current *= 0.90
    V_source[n] = V_current
    print(f"  V_source for harmonic {n}: {V_source[n]:.4f} V")
print("-" * 30)

# Step 3: Define effective capacitances for the AC model
# C_eq_CL_Cc is the equivalent series capacitance of Cc and CL
C_eq_CL_Cc = (Cc * CL) / (Cc + CL)
# C_eff is the effective capacitance for the first voltage divider stage (C_parasitic in parallel with the rest of the circuit)
C_eff = C_parasitic + C_eq_CL_Cc

# Step 4: Calculate voltages for each harmonic
V_in_peak = {}
V_L_peak = {}

print("Calculating attenuated voltages for each harmonic:")
for n in harmonics:
    f_n = n * f0
    omega_n = 2 * math.pi * f_n
    
    # Calculate frequency-dependent parasitic resistance
    R_parasitic_n = R0 * (f_n / f0)**2
    
    # --- First voltage divider: V_s -> R_p -> C_eff ---
    # Transfer function H1 = V_in / V_s
    # |H1| = 1 / sqrt(1 + (w*R*C)^2)
    H1_mag = 1 / math.sqrt(1 + (omega_n * R_parasitic_n * C_eff)**2)
    V_in_peak[n] = V_source[n] * H1_mag
    
    # --- Second voltage divider: V_in -> Cc -> CL ---
    # Transfer function H2 = V_L / V_in
    # |H2| = Z_L / (Z_Cc + Z_L) = (1/w*CL) / (1/w*Cc + 1/w*CL) = Cc / (Cc + CL)
    H2_mag = Cc / (Cc + CL)
    V_L_peak[n] = V_in_peak[n] * H2_mag
    
    print(f"Harmonic {n} (f = {f_n/1e6:.0f} MHz):")
    print(f"  R_parasitic = {R_parasitic_n:.2f} Ohms")
    print(f"  V_in peak = {V_in_peak[n]:.4f} V")
    print(f"  V_L peak (ripple) = {V_L_peak[n]:.4f} V")
print("-" * 30)

# Step 5: Calculate the DC output voltage
# Approximate V_DC as the peak of the fundamental component at the rectifier input
V_DC = V_in_peak[1]
print("Calculating DC voltage across C_L:")
print(f"V_DC is approximated by the peak input voltage of the fundamental harmonic.")
print(f"V_DC ≈ V_in_peak(n=1) = {V_DC:.4f} V")
print("-" * 30)

# Step 6: Calculate total RMS ripple voltage across C_L
sum_of_squares_ripple = sum([v**2 for v in V_L_peak.values()])
V_L_rms_total = math.sqrt(sum_of_squares_ripple / 2)
print("Calculating total RMS ripple voltage across C_L:")
print(f"V_ripple_rms^2 = (V_L1^2 + V_L3^2 + V_L5^2 + V_L7^2) / 2")
print(f"V_ripple_rms^2 = ({V_L_peak[1]:.4f}^2 + {V_L_peak[3]:.4f}^2 + {V_L_peak[5]:.4f}^2 + {V_L_peak[7]:.4f}^2) / 2")
print(f"V_ripple_rms^2 = ({sum_of_squares_ripple:.4f}) / 2 = {sum_of_squares_ripple/2:.4f}")
print(f"V_ripple_rms = {V_L_rms_total:.4f} V")
print("-" * 30)

# Step 7: Calculate the total voltage across C_L
V_total = math.sqrt(V_DC**2 + V_L_rms_total**2)
print("Calculating the total voltage across C_L:")
print(f"V_total = sqrt(V_DC^2 + V_ripple_rms^2)")
print(f"V_total = sqrt({V_DC:.4f}^2 + {V_L_rms_total:.4f}^2)")
print(f"V_total = sqrt({V_DC**2:.4f} + {V_L_rms_total**2:.4f})")
print(f"V_total = sqrt({V_DC**2 + V_L_rms_total**2:.4f})")
print(f"V_total = {V_total:.4f} V")

<<<0.9793>>>