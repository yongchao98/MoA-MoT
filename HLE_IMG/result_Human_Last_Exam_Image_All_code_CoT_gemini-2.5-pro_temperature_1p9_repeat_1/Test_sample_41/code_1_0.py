import math

# Step 1: Define constants and assumptions from the problem statement
# Given parameters
V_hvdc = 500e3  # V, HVDC voltage
L = 0.1  # H, Inductance
voltage_drop_percent = 2.5  # %
h3_distortion_percent = 5.0   # %
h5_distortion_percent = 10.0  # %
f = 50.0  # Hz, Assumed AC system frequency

# Assumptions made to solve the problem
# 1. AC line-to-line voltage (V_ac_LL) is assumed to be the same as the HVDC voltage.
V_ac_LL = V_hvdc
# 2. The given inductance L is the Thevenin equivalent inductance of the AC system at the bus.
L_th = L
# 3. The distortion percentages are for voltage harmonics (Vn / V1).

# Step 2: Calculate Reactive Power Compensation for Harmonics (Q_harmonics)
# Calculate fundamental values
omega = 2 * math.pi * f
V_ac_phase = V_ac_LL / math.sqrt(3)

# System reactances at fundamental and harmonic frequencies
X_th = omega * L_th  # Thevenin reactance at fundamental frequency (n=1)
X3 = 3 * omega * L_th # Reactance for 3rd harmonic
X5 = 5 * omega * L_th # Reactance for 5th harmonic

# Harmonic voltages
V3_phase = (h3_distortion_percent / 100) * V_ac_phase
V5_phase = (h5_distortion_percent / 100) * V_ac_phase

# Reactive power consumed by the system impedance due to harmonics (3-phase)
# Formula: Q_3ph = 3 * (V_phase^2 / X_harmonic)
Q_h3 = 3 * (V3_phase**2) / X3
Q_h5 = 3 * (V5_phase**2) / X5

# Total reactive power for harmonics
Q_harmonics = Q_h3 + Q_h5

# Step 3: Calculate Reactive Power Compensation for Voltage Drop (Q_voltage)
# The compensation is needed to bring the voltage back up by 2.5%
dV_pu = voltage_drop_percent / 100

# Formula for reactive power injection to raise voltage: dQ = (V_LL^2 / X_th) * dV_pu
Q_voltage = (V_ac_LL**2 / X_th) * dV_pu

# Step 4: Calculate Total Reactive Power Compensation (Q_total)
Q_total = Q_harmonics + Q_voltage

# Print the results in MVar
print("--- Calculation Breakdown ---")
print(f"Assumed AC Voltage (Line-to-Line): {V_ac_LL/1e3:.2f} kV")
print(f"System Thevenin Reactance (X_th): {X_th:.2f} Ohms\n")

print("1. Reactive Power for Harmonic Compensation:")
print(f"  - Compensation for 3rd Harmonic (Q_h3): {Q_h3/1e6:.2f} MVar")
print(f"  - Compensation for 5th Harmonic (Q_h5): {Q_h5/1e6:.2f} MVar")
print(f"  Total Harmonic Compensation (Q_harmonics): {Q_harmonics/1e6:.2f} MVar\n")

print("2. Reactive Power for Voltage Restoration:")
print(f"  - Compensation for {voltage_drop_percent}% Voltage Drop (Q_voltage): {Q_voltage/1e6:.2f} MVar\n")

print("3. Total Required Reactive Power Compensation:")
print(f"Total Compensation (Q_total) = Q_harmonics + Q_voltage")
print(f"Q_total = {Q_harmonics/1e6:.2f} MVar + {Q_voltage/1e6:.2f} MVar = {Q_total/1e6:.2f} MVar")

# Final answer in the required format
# <<<221.49>>>