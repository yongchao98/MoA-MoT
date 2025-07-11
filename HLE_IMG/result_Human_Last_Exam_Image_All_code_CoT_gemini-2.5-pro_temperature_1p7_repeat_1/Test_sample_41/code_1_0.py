import math

# Step 1: Define given parameters and assumptions
P_inv = 254.97  # Inverter power in MW
V_dc = 500  # HVDC voltage in kV
L = 0.1  # System inductance in H
voltage_drop_percent = 2.5 # Percentage voltage drop at Bus 6

# Harmonic distortion percentages
HD_3 = 0.05  # 5% for 3rd harmonic
HD_5 = 0.10  # 10% for 5th harmonic

# Assumptions
f = 50  # Assumed fundamental grid frequency in Hz
SCR = 3.0  # Assumed Short Circuit Ratio for the AC system

print("Problem Parameters and Assumptions:")
print(f"Inverter Power (P_inv): {P_inv} MW")
print(f"HVDC Voltage (V_dc): {V_dc} kV")
print(f"System Inductance (L): {L} H")
print(f"Voltage Drop: {voltage_drop_percent}%")
print(f"Assumed Frequency (f): {f} Hz")
print(f"Assumed Short Circuit Ratio (SCR): {SCR}")
print("-" * 30)

# Step 2: Calculate Reactive Power due to Harmonics (Q_harmonics)

# 2.1 Calculate DC current (I_d)
I_d = (P_inv * 1e6) / (V_dc * 1e3)  # in Amperes

# 2.2 Calculate fundamental AC current (I_1)
# Using the relation for a 6/12-pulse converter: I1_rms = (sqrt(6)/pi) * Id
I_1 = (math.sqrt(6) / math.pi) * I_d  # in Amperes

# 2.3 Calculate harmonic currents
I_3 = HD_3 * I_1
I_5 = HD_5 * I_1

# 2.4 Calculate harmonic reactances
omega = 2 * math.pi * f
X_3 = 3 * omega * L
X_5 = 5 * omega * L

# 2.5 Calculate reactive power consumed by harmonics (3-phase)
Q_h3 = 3 * (I_3 ** 2) * X_3 / 1e6 # in MVar
Q_h5 = 3 * (I_5 ** 2) * X_5 / 1e6 # in MVar
Q_harmonics = Q_h3 + Q_h5

print("Harmonic Compensation Calculation:")
print(f"DC Current (I_d): {I_d:.2f} A")
print(f"Fundamental AC Current (I_1): {I_1:.2f} A")
print(f"3rd Harmonic Current (I_3): {I_3:.2f} A")
print(f"5th Harmonic Current (I_5): {I_5:.2f} A")
print(f"Reactive Power for 3rd Harmonic (Q_h3): {Q_h3:.4f} MVar")
print(f"Reactive Power for 5th Harmonic (Q_h5): {Q_h5:.4f} MVar")
print(f"Total Reactive Power for Harmonics (Q_harmonics): {Q_harmonics:.4f} MVar")
print("-" * 30)

# Step 3: Calculate Reactive Power for Voltage Support (Q_fault)

# 3.1 Calculate Short Circuit Capacity (S_sc)
S_sc = SCR * P_inv  # in MVA

# 3.2 Calculate required reactive power injection
voltage_drop_pu = voltage_drop_percent / 100
Q_fault = S_sc * voltage_drop_pu  # in MVar

print("Voltage Support Calculation:")
print(f"Estimated Short Circuit Capacity (S_sc): {S_sc:.2f} MVA")
print(f"Required Reactive Power for Voltage Support (Q_fault): {Q_fault:.4f} MVar")
print("-" * 30)

# Step 4: Calculate Total Reactive Power Compensation
Q_total = Q_harmonics + Q_fault

print("Total Reactive Power Compensation Calculation:")
print("The total compensation is the sum of the compensation for harmonics and for voltage support.")
print(f"Final Equation: {Q_harmonics:.2f} MVar + {Q_fault:.2f} MVar = {Q_total:.2f} MVar")

# Final Answer
# print(f"<<<{Q_total:.2f}>>>")