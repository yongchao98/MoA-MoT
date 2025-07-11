import math

# Step 1: Define initial values from the problem description
# Load Information
P_load1 = 125  # MW
Q_load1 = 50   # MVAr

P_load2 = 90   # MW
Q_load2 = 30   # MVAr

P_load3 = 100  # MW
Q_load3 = 35   # MVAr

# Line Information
Z_line_resistance_per_km = 0.05  # Ohm/km
Z_line_reactance_per_km = 0.12 # Ohm/km
line_length = 80  # km
V_line = 230  # kV

# Support Information
S3_support_factor = 0.50  # 50% of Load 1's reactive power
Q_bess_support = 10  # MVAr

# System Factors
Q_demand_increase_factor = 0.10  # 10% increase

print("--- Step 1: Calculate Total Load Reactive Power ---")
# Sum of initial reactive loads
Q_total_initial = Q_load1 + Q_load2 + Q_load3
print(f"Initial reactive power demand from loads: {Q_load1} MVAr + {Q_load2} MVAr + {Q_load3} MVAr = {Q_total_initial:.2f} MVAr")

# Calculate increased reactive power demand
Q_load1_new = Q_load1 * (1 + Q_demand_increase_factor)
Q_load2_new = Q_load2 * (1 + Q_demand_increase_factor)
Q_load3_new = Q_load3 * (1 + Q_demand_increase_factor)
Q_total_load_new = Q_load1_new + Q_load2_new + Q_load3_new
print(f"Total reactive power demand from loads after 10% increase: {Q_total_initial:.2f} MVAr * 1.1 = {Q_total_load_new:.2f} MVAr")

print("\n--- Step 2: Calculate Reactive Power Loss in the Transmission Line ---")
# Total real power
P_total = P_load1 + P_load2 + P_load3
print(f"Total real power demand: {P_load1} MW + {P_load2} MW + {P_load3} MW = {P_total:.2f} MW")

# Total apparent power of the loads
S_total_load = math.sqrt(P_total**2 + Q_total_load_new**2)
print(f"Total apparent power (S) of loads: sqrt({P_total:.2f}^2 + {Q_total_load_new:.2f}^2) = {S_total_load:.2f} MVA")

# Calculate line current assuming total load flows through the line
# I = S / (sqrt(3) * V_L-L)
current_A = (S_total_load * 1e6) / (math.sqrt(3) * V_line * 1e3)
print(f"Calculated line current (I): {S_total_load:.2f} MVA / (sqrt(3) * {V_line} kV) = {current_A:.2f} A")

# Calculate total line reactance
X_line_total = Z_line_reactance_per_km * line_length
print(f"Total line reactance (X): {Z_line_reactance_per_km} Ohm/km * {line_length} km = {X_line_total:.2f} Ohms")

# Calculate reactive power loss in the line (Q_loss = 3 * I^2 * X)
# The result will be in VAR, so divide by 1e6 for MVAr
Q_loss_line = (3 * (current_A**2) * X_line_total) / 1e6
print(f"Reactive power loss on the line (Q_loss): 3 * {current_A:.2f}^2 A * {X_line_total:.2f} Ohms = {Q_loss_line:.2f} MVAr")

print("\n--- Step 3: Calculate Total System Reactive Demand ---")
total_q_demand = Q_total_load_new + Q_loss_line
print(f"Total Reactive Demand = (Load Demand) + (Line Loss)")
print(f"Total Reactive Demand = {Q_total_load_new:.2f} MVAr + {Q_loss_line:.2f} MVAr = {total_q_demand:.2f} MVAr")


print("\n--- Step 4: Calculate Total Available Reactive Support ---")
# Support from Generator S3
q_s3_support = Q_load1 * S3_support_factor
print(f"Reactive power support from Generator S3: {S3_support_factor*100}% of {Q_load1} MVAr = {q_s3_support:.2f} MVAr")
# Support from BESS
print(f"Reactive power support from BESS: {Q_bess_support:.2f} MVAr")

total_q_support = q_s3_support + Q_bess_support
print(f"Total available reactive support = {q_s3_support:.2f} MVAr + {Q_bess_support:.2f} MVAr = {total_q_support:.2f} MVAr")

print("\n--- Step 5: Calculate Final Reactive Power Deficit ---")
# Deficit = Total Demand - Total Support
deficit = total_q_demand - total_q_support
print(f"Reactive Power Deficit = Total Demand - Total Support")
print(f"Final Equation: ({Q_total_load_new:.2f} + {Q_loss_line:.2f}) - ({q_s3_support:.2f} + {Q_bess_support:.2f}) = {deficit:.2f}")
print(f"The total reactive power deficit is {deficit:.2f} MVAr.")
