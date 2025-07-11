import math

# Step 1: Define initial parameters from the problem description
# Load Information (P in MW, Q in MVAr)
p_load1 = 125.0
q_load1_initial = 50.0

p_load2 = 90.0
q_load2_initial = 30.0

p_load3 = 100.0
q_load3_initial = 35.0

# Line Information
z_line_reactance_per_km = 0.12  # in Ω/km
line_length_km = 80.0
line_voltage_kv = 230.0

# Voltage Drop Effect
q_demand_increase_factor = 1.10  # 10% increase

# Reactive Power Support
s3_support_factor = 0.50  # 50% of Load 1's Q
bess_support_mvar = 10.0

# Step 2: Calculate adjusted reactive power demand
q_load1_adj = q_load1_initial * q_demand_increase_factor
q_load2_adj = q_load2_initial * q_demand_increase_factor
q_load3_adj = q_load3_initial * q_demand_increase_factor
total_adj_q_demand = q_load1_adj + q_load2_adj + q_load3_adj

# Step 3: Calculate total system power for loss calculation
total_p_demand = p_load1 + p_load2 + p_load3
# Apparent power S = sqrt(P^2 + Q^2)
total_s_demand_mva_sq = total_p_demand**2 + total_adj_q_demand**2
total_s_demand_mva = math.sqrt(total_s_demand_mva_sq)

# Step 4: Calculate line reactance
total_line_reactance_ohm = z_line_reactance_per_km * line_length_km

# Step 5: Calculate reactive power loss (Q_loss = |S_MVA|^2 * X_ohm / |V_kV|^2)
q_loss_mvar = (total_s_demand_mva_sq * total_line_reactance_ohm) / (line_voltage_kv**2)

# Step 6: Calculate total reactive power requirement
total_q_requirement = total_adj_q_demand + q_loss_mvar

# Step 7: Calculate total available reactive power support
q_s3_support = q_load1_initial * s3_support_factor
total_q_support = q_s3_support + bess_support_mvar

# Step 8: Calculate the final reactive power deficit
deficit_mvar = total_q_requirement - total_q_support

# Print the results and the final equation step-by-step
print("--- Calculation of Reactive Power Deficit ---")
print("\n1. Adjusted Reactive Power Demand:")
print(f"  - Increased demand at Load 1: {q_load1_initial:.1f} MVAr * {q_demand_increase_factor-1:.0%} = {q_load1_adj:.2f} MVAr")
print(f"  - Increased demand at Load 2: {q_load2_initial:.1f} MVAr * {q_demand_increase_factor-1:.0%} = {q_load2_adj:.2f} MVAr")
print(f"  - Increased demand at Load 3: {q_load3_initial:.1f} MVAr * {q_demand_increase_factor-1:.0%} = {q_load3_adj:.2f} MVAr")
print(f"  - Total Adjusted Demand (Q_demand) = {q_load1_adj:.2f} + {q_load2_adj:.2f} + {q_load3_adj:.2f} = {total_adj_q_demand:.2f} MVAr")

print("\n2. Reactive Power Loss:")
print(f"  - Total Line Reactance (X) = {z_line_reactance_per_km:.2f} Ω/km * {line_length_km:.0f} km = {total_line_reactance_ohm:.2f} Ω")
print(f"  - Total Apparent Power (S) = sqrt({total_p_demand:.1f} MW^2 + {total_adj_q_demand:.2f} MVAr^2) = {total_s_demand_mva:.2f} MVA")
print(f"  - Reactive Power Loss (Q_loss) = ({total_s_demand_mva:.2f} MVA)^2 * {total_line_reactance_ohm:.2f} Ω / ({line_voltage_kv:.0f} kV)^2 = {q_loss_mvar:.2f} MVAr")

print("\n3. Total Reactive Power Requirement:")
print(f"  - Total Requirement = Q_demand + Q_loss = {total_adj_q_demand:.2f} MVAr + {q_loss_mvar:.2f} MVAr = {total_q_requirement:.2f} MVAr")

print("\n4. Available Reactive Power Support:")
print(f"  - Support from S3 = {s3_support_factor:.0%} of {q_load1_initial:.1f} MVAr = {q_s3_support:.2f} MVAr")
print(f"  - Support from BESS = {bess_support_mvar:.2f} MVAr")
print(f"  - Total Support = {q_s3_support:.2f} MVAr + {bess_support_mvar:.2f} MVAr = {total_q_support:.2f} MVAr")

print("\n5. Final Reactive Power Deficit Calculation:")
print(f"Deficit = (Total Requirement) - (Total Support)")
print(f"Deficit = ({total_adj_q_demand:.2f} + {q_loss_mvar:.2f}) - ({q_s3_support:.2f} + {bess_support_mvar:.2f})")
print(f"Deficit = {total_q_requirement:.2f} - {total_q_support:.2f}")
print(f"Deficit = {deficit_mvar:.2f} MVAr")

print("\nFinal Answer:")
print(f"The total reactive power deficit is {deficit_mvar:.2f} MVAr.")
<<<112.41>>>