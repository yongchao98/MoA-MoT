import cmath
import math

# Step 1: Define constants and parameters
P_in = 10e-3  # Input power in Watts (10 mW)
f = 0.8e9  # Operating frequency in Hz (0.8 GHz)
Z0 = 50.0  # Characteristic impedance in Ohms
C1 = 0.4e-12 # Capacitance of C1 in Farads (0.4 pF)
L1 = 43e-9   # Inductance of L1 in Henrys (43 nH)
C2 = 0.2e-12 # Capacitance of C2 in Farads (0.2 pF)
L2 = 39e-9   # Inductance of L2 in Henrys (39 nH)
R_L = 2.7e3  # Load resistance in Ohms (2.7 kOhm)
Q_C = 150.0  # Quality factor of capacitors
# From the graph at 800 MHz, the Quality Factor for the inductor is ~95
Q_L = 95.0

print("Step 1: System Parameters")
print(f"  Input Power (P_in): {P_in*1000} mW")
print(f"  Frequency (f): {f/1e9} GHz")
print(f"  Load Resistance (R_L): {R_L/1000} kOhm")
print(f"  Inductor Quality Factor (Q_L): {Q_L}")
print(f"  Capacitor Quality Factor (Q_C): {Q_C}")
print("-" * 30)

# Step 2: Calculate angular frequency
omega = 2 * math.pi * f

# Step 3: Calculate component loss resistances
# Loss for L1-C1 series network
X_L1 = omega * L1
R_sL1 = X_L1 / Q_L
X_C1 = 1 / (omega * C1)
R_sC1 = X_C1 / Q_C  # Note: Q_C definition gives R_s = 1 / (omega * C * Q_C), which is X_C / Q_C
R_s_L1C1 = R_sL1 + R_sC1

# Loss for L2-C2 parallel network (calculate equivalent series resistance)
X_L2 = omega * L2
R_pL2 = Q_L * X_L2
X_C2 = 1 / (omega * C2)
R_pC2 = Q_C * X_C2 # Note: Q_C definition gives R_p = Q_C / (omega * C), which is Q_C * X_C
G_L2 = 1 / R_pL2
G_C2 = 1 / R_pC2
B_L2 = -1 / X_L2
B_C2 = 1 / X_C2
Y_L2C2 = (G_L2 + G_C2) + 1j * (B_L2 + B_C2)
Z_L2C2 = 1 / Y_L2C2
R_s_L2C2_eq = Z_L2C2.real

print("Step 2: Component Loss Calculation")
print(f"  Series loss resistance of L1-C1 network: {R_s_L1C1:.3f} Ohms")
print(f"  Equivalent series loss resistance of L2-C2 network: {R_s_L2C2_eq:.3f} Ohms")
print("-" * 30)

# Step 4: Model the rectifier and calculate efficiencies
# Assumption: Ideal diode voltage doubler input resistance
R_in_diode = R_L / 2

# Total input resistance of the rectifier block
R_in_rect = R_s_L1C1 + R_s_L2C2_eq + R_in_diode

# Efficiency of the matching network
Q_match = math.sqrt(R_in_rect / Z0 - 1)
eta_match = 1 - (Q_match / Q_L)

# Efficiency of the rectifier's internal LC components
eta_LC_loss = R_in_diode / R_in_rect

# Assumption: Ideal diode RF-DC conversion efficiency
eta_diode_conv = 8 / (math.pi**2)

# Total power conversion efficiency
eta_total = eta_match * eta_LC_loss * eta_diode_conv

print("Step 3: Efficiency Calculation")
print(f"  Assumed ideal diode input resistance (R_in_diode): {R_in_diode:.1f} Ohms")
print(f"  Total rectifier input resistance (R_in_rect): {R_in_rect:.1f} Ohms")
print(f"  Matching network Q (Q_match): {Q_match:.3f}")
print(f"  Matching network efficiency (eta_match): {eta_match:.4f}")
print(f"  Rectifier LC network efficiency (eta_LC_loss): {eta_LC_loss:.4f}")
print(f"  Ideal diode conversion efficiency (eta_diode_conv): {eta_diode_conv:.4f}")
print(f"  Overall Power Conversion Efficiency (eta_total): {eta_total:.4f}")
print("-" * 30)

# Step 5: Calculate final output voltage
P_DC_out = P_in * eta_total
V_DC = math.sqrt(P_DC_out * R_L)

print("Step 4: Final Calculation")
print(f"The total efficiency of the system is estimated to be {eta_total*100:.2f}%.")
print(f"The DC power at the load is P_out = {P_in:.3f} W * {eta_total:.4f} = {P_DC_out*1000:.2f} mW.")
print(f"The final equation for the voltage across the load is:")
print(f"V_DC = sqrt(P_in * eta_total * R_L)")
print(f"V_DC = sqrt({P_in:.4f} W * {eta_total:.4f} * {R_L:.0f} Ohm) = {V_DC:.2f} V")
print(f"\nFinal Answer: {V_DC}")
