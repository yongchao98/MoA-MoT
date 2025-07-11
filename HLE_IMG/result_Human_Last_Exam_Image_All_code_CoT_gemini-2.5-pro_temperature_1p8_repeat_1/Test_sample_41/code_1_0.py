import math

# --- Given Data ---
# Power at Rectifier (Bus 5)
P_rec = 255.33  # MW
Q_rec_consumed = 74.88  # MVar

# Power from Inverter (Bus 6 to Bus 11)
P_inv = 254.97  # MW
Q_net_inv_supplied = 14.58  # MVar

# System conditions
voltage_drop_fraction = 0.025  # 2.5%
h3_distortion = 0.05  # 5%
h5_distortion = 0.10  # 10%

# --- Step 1: Calculate Reactive Power Compensation for Voltage Drop ---
print("--- Step 1: Calculate Compensation for Voltage Drop ---")

# Estimate Q/P ratio from the rectifier side
q_per_p_ratio = Q_rec_consumed / P_rec
print(f"Rectifier Q/P ratio = {Q_rec_consumed:.2f} MVar / {P_rec:.2f} MW = {q_per_p_ratio:.4f}")

# Estimate reactive power consumed by the inverter
Q_inv_consumed = q_per_p_ratio * P_inv
print(f"Estimated Inverter Q consumption = {q_per_p_ratio:.4f} * {P_inv:.2f} MW = {Q_inv_consumed:.2f} MVar")

# Estimate the nominal reactive power supplied by the capacitor bank
# It must cover the inverter's need and the net supply to the grid.
Q_cap_nominal = Q_inv_consumed + Q_net_inv_supplied
print(f"Nominal Capacitor Bank Size (Q_cap) = {Q_inv_consumed:.2f} MVar (for inverter) + {Q_net_inv_supplied:.2f} MVar (to grid) = {Q_cap_nominal:.2f} MVar")

# Calculate the loss of reactive power due to voltage drop
# Loss = Q_cap_nominal * (1 - (V_new/V_nominal)^2)
Q_volt_comp = Q_cap_nominal * (1 - (1 - voltage_drop_fraction)**2)
print(f"Compensation for Voltage Drop = Q_cap * (1 - (1 - {voltage_drop_fraction})**2)")
print(f"Q_volt_comp = {Q_cap_nominal:.2f} MVar * (1 - {1-voltage_drop_fraction}**2) = {Q_volt_comp:.2f} MVar")
print("-" * 50)


# --- Step 2: Calculate Reactive Power Compensation for Harmonics ---
print("--- Step 2: Calculate Compensation for Harmonics ---")

# Calculate fundamental apparent power (S1) at the inverter
S1 = math.sqrt(P_inv**2 + Q_net_inv_supplied**2)
print(f"Fundamental Apparent Power (S1) = sqrt({P_inv:.2f}^2 + {Q_net_inv_supplied:.2f}^2) = {S1:.2f} MVA")

# Calculate MVA rating required for 5th harmonic filter
Q_filter_5 = S1 * h5_distortion
print(f"5th Harmonic Filter Rating = S1 * {h5_distortion*100}% = {S1:.2f} MVA * {h5_distortion} = {Q_filter_5:.2f} MVar")

# Calculate MVA rating required for 3rd harmonic filter
Q_filter_3 = S1 * h3_distortion
print(f"3rd Harmonic Filter Rating = S1 * {h3_distortion*100}% = {S1:.2f} MVA * {h3_distortion} = {Q_filter_3:.2f} MVar")

# Total compensation for harmonics is the sum of filter ratings
Q_harm_comp = Q_filter_5 + Q_filter_3
print(f"Total Harmonic Compensation (Q_harm_comp) = {Q_filter_5:.2f} MVar + {Q_filter_3:.2f} MVar = {Q_harm_comp:.2f} MVar")
print("-" * 50)

# --- Step 3: Calculate Total Reactive Power Compensation ---
print("--- Step 3: Calculate Total Compensation ---")
Q_total_comp = Q_volt_comp + Q_harm_comp
print(f"Total Reactive Power Compensation = Q_volt_comp + Q_harm_comp")
print(f"Q_total = {Q_volt_comp:.2f} MVar + {Q_harm_comp:.2f} MVar = {Q_total_comp:.2f} MVar")
print("-" * 50)

final_answer = Q_total_comp
print(f"\nThe total reactive power compensation required is {final_answer:.2f} MVar.")
print(f'<<<42.72>>>')