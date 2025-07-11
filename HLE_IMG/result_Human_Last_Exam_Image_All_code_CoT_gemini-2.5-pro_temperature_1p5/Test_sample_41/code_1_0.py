import math

# --- Given Parameters ---
# From the problem description
voltage_drop_percentage = 2.5  # %
harmonic_5th_distortion = 10.0 # %
harmonic_3rd_distortion = 5.0  # %

# From the diagram at the Inverter (Bus 6), which is where compensation is needed.
P_inv = 254.97  # MW
Q_inv_1 = 14.58 # MVAr (This is the fundamental reactive power)

# --- Step 1: Calculate Reactive Power needed for Harmonic Compensation ---

# Convert percentages to decimal ratios for calculation
h5_ratio = harmonic_5th_distortion / 100.0
h3_ratio = harmonic_3rd_distortion / 100.0

# Calculate fundamental apparent power (S1) at the inverter
S1 = math.sqrt(P_inv**2 + Q_inv_1**2)

# Calculate the total apparent power (S_total) including the effect of harmonics
S_total = S1 * math.sqrt(1 + h3_ratio**2 + h5_ratio**2)

# Calculate the total reactive power (Q_total) from the total apparent power and active power
Q_total = math.sqrt(S_total**2 - P_inv**2)

# The reactive power compensation needed for harmonics is the difference
# between the total reactive power and the fundamental reactive power.
Q_harmonic = Q_total - Q_inv_1

# --- Step 2: Calculate Reactive Power needed for Voltage Support ---

# Convert voltage drop percentage to a decimal ratio
voltage_drop_ratio = voltage_drop_percentage / 100.0

# Use the approximation Q_comp = 2 * (dV/V) * S_load, where S_load is S_total.
Q_volt_comp = 2 * voltage_drop_ratio * S_total

# --- Step 3: Calculate the Total Reactive Power Compensation ---

# The total compensation is the sum of the compensation for harmonics and for the voltage drop.
Q_total_comp = Q_harmonic + Q_volt_comp

# --- Print the results with detailed steps and equations ---

print("### Calculation Breakdown ###\n")

print("--- Step 1: Reactive Power for Harmonic Compensation ---\n")
print(f"The active power at the inverter P is {P_inv:.2f} MW.")
print(f"The fundamental reactive power Q1 is {Q_inv_1:.2f} MVAr.")
print(f"Equation for fundamental apparent power: S1 = sqrt(P^2 + Q1^2)")
print(f"S1 = sqrt({P_inv:.2f}^2 + {Q_inv_1:.2f}^2) = {S1:.2f} MVA.\n")
print(f"Equation for total apparent power with harmonics: S_total = S1 * sqrt(1 + (I3/I1)^2 + (I5/I1)^2)")
print(f"S_total = {S1:.2f} * sqrt(1 + {h3_ratio:.2f}^2 + {h5_ratio:.2f}^2) = {S_total:.2f} MVA.\n")
print(f"Equation for total reactive power: Q_total = sqrt(S_total^2 - P^2)")
print(f"Q_total = sqrt({S_total:.2f}^2 - {P_inv:.2f}^2) = {Q_total:.2f} MVAr.\n")
print(f"Equation for harmonic reactive power: Q_harmonic = Q_total - Q1")
print(f"Q_harmonic = {Q_total:.2f} MVAr - {Q_inv_1:.2f} MVAr = {Q_harmonic:.2f} MVAr.\n")

print("--- Step 2: Reactive Power for Voltage Support ---\n")
print(f"The voltage drop at Bus 6 is {voltage_drop_percentage:.1f}%.")
print(f"Equation for voltage support compensation: Q_volt_comp = 2 * (ΔV/V) * S_total")
print(f"Q_volt_comp = 2 * {voltage_drop_ratio:.3f} * {S_total:.2f} MVA = {Q_volt_comp:.2f} MVAr.\n")


print("--- Step 3: Total Reactive Power Compensation ---\n")
print("Equation for total compensation: Q_total_compensation = Q_harmonic + Q_volt_comp")
print(f"Q_total_compensation = {Q_harmonic:.2f} MVAr + {Q_volt_comp:.2f} MVAr = {Q_total_comp:.2f} MVAr.\n")

print("="*40)
print(f"The total reactive power compensation required is {Q_total_comp:.2f} MVAr.")
print("="*40)
