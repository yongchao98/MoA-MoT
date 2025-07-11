import math

# --- 1. Initial Parameters ---
# Load Information (from the figure and text)
P_load1 = 125.0  # MW
Q_load1_initial = 50.0  # MVAr

P_load2 = 90.0   # MW
Q_load2_initial = 30.0   # MVAr

P_load3 = 100.0  # MW
Q_load3_initial = 35.0   # MVAr

# Line Information
line_reactance_per_km = 0.12 # Ohm/km
line_length = 80.0  # km
line_voltage_kv = 230.0 # kV

# Operational Factors
q_demand_increase_factor = 0.10 # 10%
s3_support_factor = 0.50 # 50%
bess_support_mvar = 10.0 # MVAr

# --- 2. Calculations ---

# Step 2a: Calculate final reactive power demand
Q_load1_final = Q_load1_initial * (1 + q_demand_increase_factor)
Q_load2_final = Q_load2_initial * (1 + q_demand_increase_factor)
Q_load3_final = Q_load3_initial * (1 + q_demand_increase_factor)
total_q_demand = Q_load1_final + Q_load2_final + Q_load3_final

# Step 2b: Calculate reactive power loss
total_p_demand = P_load1 + P_load2 + P_load3
total_s_demand_mva = math.sqrt(total_p_demand**2 + total_q_demand**2)
total_line_reactance_ohm = line_reactance_per_km * line_length
q_loss_mvar = (total_s_demand_mva**2 * total_line_reactance_ohm) / (line_voltage_kv**2)

# Step 2c: Calculate total reactive power required
total_q_required = total_q_demand + q_loss_mvar

# Step 2d: Calculate specified reactive power supply
q_s3_supply = s3_support_factor * Q_load1_initial
total_q_supply = q_s3_supply + bess_support_mvar

# Step 2e: Calculate reactive power deficit
deficit_mvar = total_q_required - total_q_supply

# --- 3. Output ---

print("Calculating the Reactive Power Deficit")
print("=======================================\n")

# Part 1: Total Reactive Power Demand
print("--- Part 1: Total Reactive Power Demand ---")
print(f"Initial demands: Load 1 = {Q_load1_initial:.1f} MVAr, Load 2 = {Q_load2_initial:.1f} MVAr, Load 3 = {Q_load3_initial:.1f} MVAr.")
print(f"A {q_demand_increase_factor*100:.0f}% increase occurs due to voltage drop.")
print(f"Final Reactive Power Demand = ({Q_load1_initial:.1f} * 1.1) + ({Q_load2_initial:.1f} * 1.1) + ({Q_load3_initial:.1f} * 1.1)")
print(f"                          = {Q_load1_final:.1f} MVAr + {Q_load2_final:.1f} MVAr + {Q_load3_final:.2f} MVAr = {total_q_demand:.2f} MVAr\n")

# Part 2: Reactive Power Losses
print("--- Part 2: Reactive Power Losses ---")
print(f"To find line losses, we first calculate the total power flow.")
print(f"Total Real Power P = {P_load1:.1f} + {P_load2:.1f} + {P_load3:.1f} = {total_p_demand:.1f} MW")
print(f"Total Reactive Power Q = {total_q_demand:.2f} MVAr")
print(f"Total Apparent Power S = sqrt({total_p_demand:.1f}^2 + {total_q_demand:.2f}^2) = {total_s_demand_mva:.2f} MVA")
print(f"Total Line Reactance X = {line_reactance_per_km:.2f} Ohm/km * {line_length:.1f} km = {total_line_reactance_ohm:.2f} Ohms")
print(f"Reactive Power Loss Q_loss = S^2 * X / V^2 = ({total_s_demand_mva:.2f} MVA)^2 * {total_line_reactance_ohm:.2f} Ohms / ({line_voltage_kv:.1f} kV)^2 = {q_loss_mvar:.2f} MVAr\n")

# Part 3: Total Specified Supply
print("--- Part 3: Specified Reactive Power Supply ---")
print(f"Supply from S3 = {s3_support_factor*100:.0f}% of initial Q_Load1 = {s3_support_factor:.2f} * {Q_load1_initial:.1f} = {q_s3_supply:.2f} MVAr")
print(f"Supply from BESS = {bess_support_mvar:.1f} MVAr")
print(f"Total Specified Supply = {q_s3_supply:.2f} MVAr + {bess_support_mvar:.1f} MVAr = {total_q_supply:.2f} MVAr\n")

# Part 4: Final Deficit Calculation
print("--- Part 4: Final Reactive Power Deficit Calculation ---")
print(f"Total Reactive Power Required = {total_q_demand:.2f} MVAr (Demand) + {q_loss_mvar:.2f} MVAr (Loss) = {total_q_required:.2f} MVAr")
print(f"Total Specified Reactive Power Supplied = {total_q_supply:.2f} MVAr")
print(f"\nFinal Deficit = Total Required - Total Supplied")
print(f"Final Deficit = (Total Demand + Loss) - (Total Supply)")
print(f"Final Deficit = (({Q_load1_final:.1f} + {Q_load2_final:.1f} + {Q_load3_final:.2f}) + {q_loss_mvar:.2f}) - ({q_s3_supply:.2f} + {bess_support_mvar:.1f})")
print(f"Final Deficit = ({total_q_demand:.2f} + {q_loss_mvar:.2f}) - {total_q_supply:.2f}")
print(f"Final Deficit = {total_q_required:.2f} - {total_q_supply:.2f} = {deficit_mvar:.2f} MVAr")
