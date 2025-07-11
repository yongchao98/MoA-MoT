import math

# --- Given Parameters and Assumptions ---

# Active power at the inverter from the diagram
P_inv = 254.97  # MW

# Voltage drop at Bus 6 after the fault
voltage_drop_percent = 2.5  # %

# Harmonic distortion percentages
fifth_harmonic_distortion_percent = 10.0  # %
third_harmonic_distortion_percent = 5.0  # %

# System inductance and HVDC voltage
L = 0.1  # H
V_dc = 500.0  # kV

# --- Assumptions for Calculation ---
f = 50.0  # Hz (Standard system frequency)
V_ac_ll = 400.0  # kV (Assumed nominal AC line-to-line voltage at inverter bus)

# --- Step 1: Calculate Reactive Power for Voltage Support (Q_voltage) ---

# Convert percentages to decimal form
voltage_drop_frac = voltage_drop_percent / 100.0

# Calculate angular frequency
omega = 2 * math.pi * f

# Calculate the equivalent Thevenin reactance of the AC system
X_sys = omega * L

# Calculate the nominal AC phase voltage
V_ac_ph_kv = V_ac_ll / math.sqrt(3)
V_ac_ph_v = V_ac_ph_kv * 1000

# The formula for reactive power injection to cause a voltage change is dQ = dV * V / X
# For a desired voltage change of 'voltage_drop_frac', the required reactive power per phase is:
# Q_phase = (voltage_drop_frac * V_ac_ph) * V_ac_ph / X_sys
Q_voltage_per_phase_Mvar = (voltage_drop_frac * V_ac_ph_v**2) / X_sys / 1e6

# Total three-phase reactive power for voltage support
Q_voltage = 3 * Q_voltage_per_phase_Mvar

print("--- Calculation for Reactive Power Compensation ---")
print("\nStep 1: Calculate compensation for 2.5% voltage drop")
print(f"Assumed AC Voltage (Line-to-Line): {V_ac_ll:.2f} kV")
print(f"System Reactance (X = 2*pi*f*L): {X_sys:.2f} Ohms")
print(f"Required compensation for voltage support (Q_voltage): {Q_voltage:.2f} MVar")

# --- Step 2: Calculate Reactive Power for Harmonic Distortion (Q_harmonics) ---

# Convert distortion percentages to decimal form
fifth_harmonic_frac = fifth_harmonic_distortion_percent / 100.0
third_harmonic_frac = third_harmonic_distortion_percent / 100.0

# Calculate reactive power based on harmonic distortion percentages of active power
Q_h5 = fifth_harmonic_frac * P_inv
Q_h3 = third_harmonic_frac * P_inv

# Total reactive power for harmonics
Q_harmonics = Q_h5 + Q_h3

print("\nStep 2: Calculate compensation for harmonic distortion")
print(f"Inverter Active Power (P_inv): {P_inv:.2f} MW")
print(f"Reactive power for 5th harmonic (10% of P_inv): {Q_h5:.2f} MVar")
print(f"Reactive power for 3rd harmonic (5% of P_inv): {Q_h3:.2f} MVar")
print(f"Total compensation for harmonics (Q_harmonics): {Q_harmonics:.2f} MVar")

# --- Step 3: Calculate Total Reactive Power Compensation ---

# Total compensation is the sum of the two components
Q_total = Q_voltage + Q_harmonics

print("\nStep 3: Calculate Total Required Compensation")
print("\nThe final equation for total reactive power compensation is:")
print(f"Q_total = Q_voltage + Q_harmonics")
print(f"Q_total = {Q_voltage:.2f} MVar + {Q_harmonics:.2f} MVar")
print(f"Total Reactive Power Compensation = {Q_total:.2f} MVar")

# Final Answer
final_answer = Q_total