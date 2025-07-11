import math

# --- Initial Parameters ---
# Load Information (initial)
P_L1 = 125.0  # MW
Q_L1_initial = 50.0  # MVAr

P_L2 = 90.0   # MW
Q_L2_initial = 30.0  # MVAr

P_L3 = 100.0  # MW
Q_L3_initial = 35.0  # MVAr

# Line Information
line_length_total = 80.0  # km
line_reactance_per_km = 0.12  # Ohm/km
V_line = 230.0  # kV

# System Adjustments and Support
q_increase_factor = 1.10  # 10% increase
s3_support_factor = 0.50   # 50% of Load 1's Q
bess_support = 10.0      # MVAr

# --- Step 1: Calculate Final Reactive Power Demand ---
print("--- Step 1: Calculating Final Reactive Power Demand ---")
Q_L1_final = Q_L1_initial * q_increase_factor
Q_L2_final = Q_L2_initial * q_increase_factor
Q_L3_final = Q_L3_initial * q_increase_factor
Q_demand_total = Q_L1_final + Q_L2_final + Q_L3_final

print(f"Final Reactive Demand for Load 1: {Q_L1_initial:.2f} MVAr * {q_increase_factor:.2f} = {Q_L1_final:.2f} MVAr")
print(f"Final Reactive Demand for Load 2: {Q_L2_initial:.2f} MVAr * {q_increase_factor:.2f} = {Q_L2_final:.2f} MVAr")
print(f"Final Reactive Demand for Load 3: {Q_L3_initial:.2f} MVAr * {q_increase_factor:.2f} = {Q_L3_final:.2f} MVAr")
print(f"Total Final Reactive Power Demand = {Q_L1_final:.2f} + {Q_L2_final:.2f} + {Q_L3_final:.2f} = {Q_demand_total:.2f} MVAr\n")


# --- Step 2: Calculate Reactive Power Line Loss ---
print("--- Step 2: Calculating Reactive Power Line Loss ---")
# Assumption: Bus 8 is halfway, so T5 and T6 are each 40 km.
line_length_segment = line_length_total / 2.0
X_segment = line_reactance_per_km * line_length_segment

# Loss in segment T6 (Bus 8 to 9), serving Load 2
S_L2_final_sq = P_L2**2 + Q_L2_final**2
Q_loss_T6 = X_segment * S_L2_final_sq / (V_line**2)

print(f"Reactance of one line segment (40km): {X_segment:.2f} Ohms")
print(f"Apparent power squared for Load 2 (S_L2^2): {P_L2**2:.2f} + {Q_L2_final**2:.2f} = {S_L2_final_sq:.2f} (MVA)^2")
print(f"Reactive Loss in T6 = {X_segment:.2f} * {S_L2_final_sq:.2f} / {V_line**2:.2f} = {Q_loss_T6:.2f} MVAr")

# Loss in segment T5 (Bus 7 to 8), serving Load 2 and Load 3
P_T5 = P_L2 + P_L3
Q_T5 = Q_L2_final + Q_L3_final
S_T5_final_sq = P_T5**2 + Q_T5**2
Q_loss_T5 = X_segment * S_T5_final_sq / (V_line**2)

print(f"Apparent power squared for T5 flow (S_T5^2): {P_T5**2:.2f} + {Q_T5**2:.2f} = {S_T5_final_sq:.2f} (MVA)^2")
print(f"Reactive Loss in T5 = {X_segment:.2f} * {S_T5_final_sq:.2f} / {V_line**2:.2f} = {Q_loss_T5:.2f} MVAr")

Q_loss_total = Q_loss_T5 + Q_loss_T6
print(f"Total Reactive Power Line Loss = {Q_loss_T5:.2f} + {Q_loss_T6:.2f} = {Q_loss_total:.2f} MVAr\n")


# --- Step 3: Calculate Total Reactive Power Required ---
print("--- Step 3: Calculating Total Reactive Power Required ---")
Q_required_total = Q_demand_total + Q_loss_total
print(f"Total Q Required = Total Q Demand + Total Q Loss = {Q_demand_total:.2f} MVAr + {Q_loss_total:.2f} MVAr = {Q_required_total:.2f} MVAr\n")


# --- Step 4: Calculate Available Reactive Power Supply ---
print("--- Step 4: Calculating Available Reactive Power Supply ---")
Q_supply_S3 = Q_L1_final * s3_support_factor
Q_supply_BESS = bess_support
Q_supply_total = Q_supply_S3 + Q_supply_BESS

print(f"Q Supply from S3 (50% of Load 1's final Q): {Q_L1_final:.2f} MVAr * {s3_support_factor:.2f} = {Q_supply_S3:.2f} MVAr")
print(f"Q Supply from BESS: {Q_supply_BESS:.2f} MVAr")
print(f"Total Available Q Supply = {Q_supply_S3:.2f} MVAr + {Q_supply_BESS:.2f} MVAr = {Q_supply_total:.2f} MVAr\n")


# --- Step 5: Calculate Reactive Power Deficit ---
print("--- Step 5: Calculating Final Reactive Power Deficit ---")
deficit = Q_required_total - Q_supply_total
print("Deficit = Total Required Q - Total Available Q")
print(f"Deficit = ({Q_demand_total:.2f} MVAr [Loads] + {Q_loss_total:.2f} MVAr [Loss]) - ({Q_supply_S3:.2f} MVAr [S3] + {Q_supply_BESS:.2f} MVAr [BESS])")
print(f"Deficit = {Q_required_total:.2f} MVAr - {Q_supply_total:.2f} MVAr = {deficit:.2f} MVAr")

print("\n<<<93.57>>>")