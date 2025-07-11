import math

# --- Given Parameters from the problem and diagram ---
# Power at the inverter (Bus 6)
P_inverter_MW = 254.97  # Active Power in MW
Q_inverter_net_MVar = 14.58 # Net Reactive Power supplied from Bus 6 in MVar

# System disturbances
voltage_drop_percent = 2.5
harmonic_distortion_5th_percent = 10
harmonic_distortion_3rd_percent = 5

# --- Step 1: Calculate compensation for harmonic distortion ---
# The required reactive power is estimated as a percentage of the active power.
# The total percentage is the sum of the given harmonic distortions.
total_harmonic_distortion_percent = harmonic_distortion_5th_percent + harmonic_distortion_3rd_percent
Q_comp_harmonics_MVar = (total_harmonic_distortion_percent / 100) * P_inverter_MW

print("--- Calculation for Required Reactive Power Compensation ---")
print("\nStep 1: Calculate compensation for harmonic distortion (Q_harm).")
print(f"This is estimated as the total harmonic distortion percentage applied to the inverter's active power.")
print(f"Q_harm = ({harmonic_distortion_5th_percent}% + {harmonic_distortion_3rd_percent}%) * {P_inverter_MW} MW")
print(f"Q_harm = {total_harmonic_distortion_percent / 100} * {P_inverter_MW} MW = {Q_comp_harmonics_MVar:.2f} MVar")

# --- Step 2: Calculate compensation for voltage support ---
# First, calculate the apparent power (S) at the inverter bus.
S_inverter_MVA = math.sqrt(P_inverter_MW**2 + Q_inverter_net_MVar**2)

# The required reactive power for voltage support is estimated as a percentage of the apparent power.
Q_comp_voltage_MVar = (voltage_drop_percent / 100) * S_inverter_MVA

print("\nStep 2: Calculate compensation for voltage support (Q_volt).")
print(f"First, find the apparent power at the inverter (S_inv).")
print(f"S_inv = sqrt({P_inverter_MW}^2 + {Q_inverter_net_MVar}^2) = {S_inverter_MVA:.2f} MVA")
print(f"The compensation is the voltage drop percentage applied to the apparent power.")
print(f"Q_volt = {voltage_drop_percent}% * {S_inverter_MVA:.2f} MVA")
print(f"Q_volt = {voltage_drop_percent / 100} * {S_inverter_MVA:.2f} MVA = {Q_comp_voltage_MVar:.2f} MVar")

# --- Step 3: Calculate Total Reactive Power Compensation ---
Q_total_MVar = Q_comp_harmonics_MVar + Q_comp_voltage_MVar

print("\nStep 3: Calculate the total reactive power compensation (Q_total).")
print("Q_total = Q_harm + Q_volt")
print(f"Q_total = {Q_comp_harmonics_MVar:.2f} MVar + {Q_comp_voltage_MVar:.2f} MVar = {Q_total_MVar:.2f} MVar")

print("\n--- Final Answer ---")
print(f"The total reactive power compensation required is {Q_total_MVar:.2f} MVar.")
<<<44.63>>>