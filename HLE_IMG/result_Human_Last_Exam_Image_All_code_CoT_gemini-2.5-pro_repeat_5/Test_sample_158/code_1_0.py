import cmath
import math

# Step 1: Define constants and parameters
f = 0.8e9  # Operating frequency in Hz
P_in = 10e-3  # Input power in W
C1_val = 0.4e-12  # C1 in F
C2_val = 0.2e-12  # C2 in F
L1_val = 43e-9  # L1 in H
L2_val = 39e-9  # L2 in H
RL = 2.7e3  # Load resistance in Ohms
QC = 150  # Quality factor of capacitors
# From the graph at 800 MHz, the Quality factor of the inductor is ~95
QL = 95

# Step 2: Calculate angular frequency and model the rectifier
omega = 2 * math.pi * f
# Model the rectifier input resistance as RL/2 for a shunt topology
R_eff = RL / 2

# Step 3: Calculate loss resistances and reactances
Rs_L1 = (omega * L1_val) / QL
X_L1 = omega * L1_val
Rs_C1 = 1 / (omega * C1_val * QC)
X_C1 = -1 / (omega * C1_val)

Rs_L2 = (omega * L2_val) / QL
X_L2 = omega * L2_val
Rs_C2 = 1 / (omega * C2_val * QC)
X_C2 = -1 / (omega * C2_val)

# Step 4: Calculate branch impedances
# Series branch (L1, C1, R_eff)
Z_s = (Rs_L1 + Rs_C1 + R_eff) + 1j * (X_L1 + X_C1)
Y_s = 1 / Z_s

# Parallel branch (L2, C2)
Z_L2 = Rs_L2 + 1j * X_L2
Z_C2 = Rs_C2 + 1j * X_C2
Y_p = 1 / Z_L2 + 1 / Z_C2

# Step 5: Calculate power distribution
# Total admittance of the rectifier circuitry
Y_tot = Y_s + Y_p
# Total power is Pin = V_rms^2 * Re(Y_tot)
# V_rms^2 is the square of the RMS voltage at the rectifier circuit input
V_rms_sq = P_in / Y_tot.real

# Power delivered to the series branch
P_s = V_rms_sq * Y_s.real

# Step 6: Calculate RF power delivered to the rectifier
# Total series resistance of the series branch
R_s_total = Z_s.real
# Power delivered to R_eff (the rectifier)
P_rect_RF = P_s * (R_eff / R_s_total)

# Step 7: Assume ideal RF to DC conversion
P_DC = P_rect_RF

# Step 8: Calculate the final DC voltage
V_DC = math.sqrt(P_DC * RL)

# Print the results and the final equation
print(f"--- Calculation Steps ---")
print(f"Operating Frequency (f): {f/1e9:.2f} GHz")
print(f"Angular Frequency (w): {omega:.2e} rad/s")
print(f"Input Power (Pin): {P_in * 1000:.1f} mW")
print(f"Load Resistor (RL): {RL/1000:.1f} kOhms")
print(f"Inductor Q (QL): {QL}")
print(f"Capacitor Q (QC): {QC}")
print("\n--- Rectifier Model ---")
print(f"Effective RF resistance of rectifier (Reff = RL/2): {R_eff:.1f} Ohms")
print("\n--- Passive Component Losses ---")
print(f"Loss resistance of L1 (Rs_L1): {Rs_L1:.3f} Ohms")
print(f"Loss resistance of C1 (Rs_C1): {Rs_C1:.3f} Ohms")
print(f"Loss resistance of L2 (Rs_L2): {Rs_L2:.3f} Ohms")
print(f"Loss resistance of C2 (Rs_C2): {Rs_C2:.3f} Ohms")
print("\n--- Power Analysis ---")
print(f"Total power delivered to rectifier circuit: {P_in * 1000:.2f} mW")
print(f"RF power delivered to rectifier (P_rect_RF): {P_rect_RF * 1000:.2f} mW")
print(f"DC Power at load (assuming ideal conversion, P_DC = P_rect_RF): {P_DC * 1000:.2f} mW")
print("\n--- Final Calculation ---")
print(f"Voltage across load (V_DC) = sqrt(P_DC * RL)")
print(f"V_DC = sqrt({P_DC:.4f} W * {RL:.1f} Ohms)")
print(f"V_DC = {V_DC:.2f} V")

print(f"\nFinal Answer:")
print(f"The voltage across the load R_L is calculated as V_DC = sqrt(P_DC * RL) = sqrt({P_DC:.4f} * {RL:.0f}) = {V_DC:.2f} V.")
<<<4.98>>>