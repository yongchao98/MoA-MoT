import math

# --- Given Parameters ---
# Real power at the inverter (Bus 6) in MW
P_inv = 254.97

# Voltage drop at Bus 6 after the fault
voltage_drop_pct = 2.5

# Harmonic distortion percentages
hd5_pct = 10.0
hd3_pct = 5.0

# --- Step-by-Step Calculation ---

# 1. Calculate the reactive power compensation required for harmonic distortion.
# This is assumed to be a percentage of the real power.
total_hd_pct = hd5_pct + hd3_pct
q_harmonics = (total_hd_pct / 100) * P_inv

# 2. Calculate the reactive power compensation required to restore voltage stability.
# A voltage drop implies a reactive power deficit. A common approximation is that the
# required reactive power injection is proportional to the voltage drop percentage
# multiplied by the real power flow.
q_voltage = (voltage_drop_pct / 100) * P_inv

# 3. Calculate the total reactive power compensation required.
# This is the sum of the compensation for harmonics and for voltage stability.
q_total = q_harmonics + q_voltage

# --- Output the results ---

print("Step 1: Calculate Reactive Power Compensation for Harmonics (Q_harm)")
print(f"Total Harmonic Distortion = {hd5_pct}% (5th) + {hd3_pct}% (3rd) = {total_hd_pct:.2f}%")
print(f"Real Power at Inverter (P_inv) = {P_inv} MW")
print(f"Q_harm = ({total_hd_pct:.2f} / 100) * {P_inv} MW = {q_harmonics:.2f} MVar\n")

print("Step 2: Calculate Reactive Power Compensation for Voltage Stability (Q_volt)")
print(f"Voltage Drop = {voltage_drop_pct:.2f}%")
print(f"Q_volt = ({voltage_drop_pct:.2f} / 100) * {P_inv} MW = {q_voltage:.2f} MVar\n")

print("Step 3: Calculate Total Reactive Power Compensation (Q_total)")
print("Q_total = Q_harm + Q_volt")
print(f"Q_total = {q_harmonics:.2f} MVar + {q_voltage:.2f} MVar = {q_total:.2f} MVar\n")

print("--- Final Equation ---")
print(f"Total Compensation = (({hd5_pct} / 100) + ({hd3_pct} / 100)) * {P_inv} + ({voltage_drop_pct} / 100) * {P_inv}")
print(f"Total Compensation = ({total_hd_pct/100}) * {P_inv} + ({voltage_drop_pct/100}) * {P_inv}")
print(f"Total Compensation = {q_harmonics:.2f} + {q_voltage:.2f} = {q_total:.2f} MVar")

print(f"\nThe total reactive power compensation required is {q_total:.2f} MVar.")
<<<44.62>>>