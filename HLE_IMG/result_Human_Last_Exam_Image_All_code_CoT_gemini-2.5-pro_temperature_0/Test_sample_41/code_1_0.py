import math

# Step 1: Define the given parameters from the problem description and diagram.
P_inv = 254.97  # Active power from inverter at Bus 6 (MW)
Q_net = 14.58   # Net reactive power injected into the grid from Bus 6 (MVar)
volt_drop_pct = 2.5  # Voltage drop at Bus 6 in percent
h3_dist_pct = 5.0    # 3rd harmonic distortion in percent
h5_dist_pct = 10.0   # 5th harmonic distortion in percent

# Assumed factor for reactive power consumption by the LCC inverter
k_inv_consumption_factor = 0.55

# Step 2: Calculate the reactive power compensation required due to the voltage drop.
# Estimate the reactive power consumed by the inverter.
Q_inv_consumed = P_inv * k_inv_consumption_factor

# Estimate the initial reactive power supplied by the capacitor bank at nominal voltage.
Q_cap_initial = Q_inv_consumed + Q_net

# Calculate the voltage ratio after the drop.
voltage_ratio = 1 - (volt_drop_pct / 100)

# Calculate the new reactive power supplied by the capacitor bank at the reduced voltage.
Q_cap_final = Q_cap_initial * (voltage_ratio ** 2)

# The compensation needed is the loss in reactive power from the capacitor bank.
Q_comp_volt = Q_cap_initial - Q_cap_final

print("--- Calculation for Voltage Drop Compensation ---")
print(f"Active Power (P): {P_inv} MW")
print(f"Net Reactive Power Output (Q_net): {Q_net} MVar")
print(f"Estimated Inverter Reactive Power Consumption: {P_inv:.2f} MW * {k_inv_consumption_factor} = {Q_inv_consumed:.2f} MVar")
print(f"Estimated Initial Capacitor Bank Size: {Q_inv_consumed:.2f} MVar + {Q_net:.2f} MVar = {Q_cap_initial:.2f} MVar")
print(f"Capacitor Output after {volt_drop_pct}% voltage drop: {Q_cap_initial:.2f} MVar * ({voltage_ratio})^2 = {Q_cap_final:.2f} MVar")
print(f"Reactive Power Compensation for Voltage Drop (Q_comp_volt): {Q_cap_initial:.2f} MVar - {Q_cap_final:.2f} MVar = {Q_comp_volt:.2f} MVar\n")


# Step 3: Calculate the reactive power compensation required due to harmonic distortion.
# Calculate the fundamental apparent power (S1).
S_fundamental = math.sqrt(P_inv**2 + Q_net**2)

# Convert harmonic distortion percentages to decimal form.
h3_dist = h3_dist_pct / 100
h5_dist = h5_dist_pct / 100

# Calculate the Total Harmonic Distortion (THD) of the current.
THD = math.sqrt(h3_dist**2 + h5_dist**2)

# Calculate the distortion power (D), which represents the reactive power due to harmonics.
Q_comp_harmonic = S_fundamental * THD

print("--- Calculation for Harmonic Distortion Compensation ---")
print(f"Fundamental Apparent Power (S1): sqrt({P_inv}^2 + {Q_net}^2) = {S_fundamental:.2f} MVA")
print(f"Total Harmonic Distortion (THD): sqrt({h3_dist}^2 + {h5_dist}^2) = {THD:.4f}")
print(f"Reactive Power Compensation for Harmonics (Q_comp_harmonic): {S_fundamental:.2f} MVA * {THD:.4f} = {Q_comp_harmonic:.2f} MVar\n")


# Step 4: Calculate the total reactive power compensation.
Q_total = Q_comp_volt + Q_comp_harmonic

print("--- Total Reactive Power Compensation ---")
print("The total reactive power compensation is the sum of the compensation for the voltage drop and for the harmonic distortion.")
print(f"Total Compensation = Q_comp_volt + Q_comp_harmonic")
print(f"Total Compensation = {Q_comp_volt:.2f} MVar + {Q_comp_harmonic:.2f} MVar = {Q_total:.2f} MVar")

# Final Answer
# print(f"\n<<<The final answer is {Q_total:.2f}>>>")