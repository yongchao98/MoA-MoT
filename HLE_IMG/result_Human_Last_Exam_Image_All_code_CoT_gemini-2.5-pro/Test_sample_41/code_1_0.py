import math

# --- Given Data and Assumptions ---

# Active power transferred by the inverter (from diagram)
P_inv = 254.97  # MW

# Reactive power supplied by the AC grid to the inverter station (from diagram)
Q_grid = 14.58  # MVAr

# Voltage drop at Bus 6 after the fault
voltage_drop_percent = 2.5  # %

# Harmonic distortion percentages given in the problem
harmonic_3rd_percent = 5.0   # %
harmonic_5th_percent = 10.0  # %

# Assumption: Reactive power consumed by an LCC converter is ~50% of active power
converter_q_factor = 0.5

# --- Calculation ---

print("Calculating the total reactive power compensation required.\n")

# Step 1: Calculate the reactive power lost due to the voltage drop
print("--- Part 1: Compensation for Voltage Drop ---")

# Estimate the total reactive power consumed by the inverter
q_converter_demand = P_inv * converter_q_factor

# Determine the initial reactive power supplied by the filter capacitor bank (C11)
q_filter_initial = q_converter_demand - Q_grid
print(f"Estimated reactive power demand of the inverter: {P_inv:.2f} MW * {converter_q_factor} = {q_converter_demand:.2f} MVAr")
print(f"Initial reactive power from filter bank (C11) = {q_converter_demand:.2f} MVAr - {Q_grid:.2f} MVAr = {q_filter_initial:.2f} MVAr\n")

# Calculate the reduction in reactive power supply from the filter due to the voltage drop
# Q is proportional to V^2. New voltage V_new = V_initial * (1 - 0.025)
v_factor = 1 - (voltage_drop_percent / 100.0)
q_filter_after_fault = q_filter_initial * (v_factor ** 2)
compensation_for_voltage_drop = q_filter_initial - q_filter_after_fault

print(f"Voltage drops by {voltage_drop_percent}%, so the new voltage is {v_factor * 100} % of the original.")
print(f"Reactive power supplied by the filter after the fault = {q_filter_initial:.2f} MVAr * ({v_factor:.3f})^2 = {q_filter_after_fault:.2f} MVAr")
print(f"Compensation needed for voltage drop = {q_filter_initial:.2f} MVAr - {q_filter_after_fault:.2f} MVAr = {compensation_for_voltage_drop:.2f} MVAr\n")


# Step 2: Calculate the reactive power required for harmonic compensation
print("--- Part 2: Compensation for Harmonic Distortion ---")

# Interpret harmonic distortion as a reactive power requirement based on active power
compensation_for_3rd_harmonic = (harmonic_3rd_percent / 100.0) * P_inv
compensation_for_5th_harmonic = (harmonic_5th_percent / 100.0) * P_inv
total_harmonic_compensation = compensation_for_3rd_harmonic + compensation_for_5th_harmonic

print(f"Compensation for 3rd harmonic ({harmonic_3rd_percent}%) = {harmonic_3rd_percent/100.0:.2f} * {P_inv:.2f} MW = {compensation_for_3rd_harmonic:.2f} MVAr")
print(f"Compensation for 5th harmonic ({harmonic_5th_percent}%) = {harmonic_5th_percent/100.0:.2f} * {P_inv:.2f} MW = {compensation_for_5th_harmonic:.2f} MVAr")
print(f"Total compensation for harmonics = {compensation_for_3rd_harmonic:.2f} MVAr + {compensation_for_5th_harmonic:.2f} MVAr = {total_harmonic_compensation:.2f} MVAr\n")


# Step 3: Calculate the total reactive power compensation
print("--- Part 3: Total Compensation ---")

total_compensation = compensation_for_voltage_drop + total_harmonic_compensation

print("The total reactive power compensation required is the sum of the two components.")
print(f"Total Compensation = (Compensation for Voltage Drop) + (Compensation for Harmonics)")
print(f"Total Compensation = {compensation_for_voltage_drop:.2f} MVAr + {total_harmonic_compensation:.2f} MVAr")
print(f"Total Compensation = {total_compensation:.2f} MVAr")
print("\n<<<43.82>>>")
