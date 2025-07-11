import math

# 1. Define initial parameters from the problem description
# Load Information
P_load1 = 125.0  # MW
Q_load1_initial = 50.0  # MVAr

P_load2 = 90.0  # MW
Q_load2_initial = 30.0  # MVAr

P_load3 = 100.0  # MW
Q_load3_initial = 35.0  # MVAr

# Line Information
line_impedance_x_per_km = 0.12  # Ohm/km
line_length = 80.0  # km
line_voltage = 230.0  # kV

# Factors
q_demand_increase_factor = 0.10  # 10% increase in reactive power demand

# Reactive Power Support
s3_support_factor = 0.50  # 50% of Load 1's reactive power
battery_support = 10.0  # MVAr

# --- Calculations ---

# Step 1: Calculate the total adjusted reactive power demand from loads
print("### Step 1: Calculate Adjusted Reactive Power Demand ###")
q_demand_initial = Q_load1_initial + Q_load2_initial + Q_load3_initial
q_demand_adjusted = q_demand_initial * (1 + q_demand_increase_factor)
print(f"Total initial reactive demand from loads = {Q_load1_initial} + {Q_load2_initial} + {Q_load3_initial} = {q_demand_initial:.2f} MVAr")
print(f"Adjusted reactive demand = {q_demand_initial:.2f} * (1 + {q_demand_increase_factor}) = {q_demand_adjusted:.2f} MVAr\n")

# Step 2: Calculate the reactive power loss in the transmission line
print("### Step 2: Calculate Reactive Power Loss ###")
p_total = P_load1 + P_load2 + P_load3
s_total_sq = p_total**2 + q_demand_adjusted**2
x_line_total = line_impedance_x_per_km * line_length
q_loss_line = (s_total_sq * x_line_total) / (line_voltage**2)
print(f"Total real power P = {P_load1} + {P_load2} + {P_load3} = {p_total:.2f} MW")
print(f"Total apparent power squared S^2 = {p_total:.2f}^2 + {q_demand_adjusted:.2f}^2 = {s_total_sq:.2f} MVA^2")
print(f"Total line reactance X = {line_impedance_x_per_km} Ω/km * {line_length} km = {x_line_total:.2f} Ω")
print(f"Reactive Power Loss Q_loss = S^2 * X / V^2 = {s_total_sq:.2f} * {x_line_total:.2f} / {line_voltage**2} = {q_loss_line:.2f} MVAr\n")

# Step 3: Calculate the total reactive power need
print("### Step 3: Calculate Total Reactive Power Need ###")
q_total_need = q_demand_adjusted + q_loss_line
print(f"Total Reactive Need = Adjusted Demand + Line Loss")
print(f"Total Reactive Need = {q_demand_adjusted:.2f} + {q_loss_line:.2f} = {q_total_need:.2f} MVAr\n")

# Step 4: Calculate the total available reactive power supply
print("### Step 4: Calculate Total Available Reactive Power Supply ###")
q_supply_s3 = Q_load1_initial * s3_support_factor
q_total_supply = q_supply_s3 + battery_support
print(f"Supply from S3 = {s3_support_factor:.2f} * {Q_load1_initial} = {q_supply_s3:.2f} MVAr")
print(f"Supply from Battery = {battery_support:.2f} MVAr")
print(f"Total Available Supply = {q_supply_s3:.2f} + {battery_support:.2f} = {q_total_supply:.2f} MVAr\n")

# Step 5: Calculate the final reactive power deficit
print("### Step 5: Calculate Final Reactive Power Deficit ###")
deficit = q_total_need - q_total_supply
print("Deficit = Total Need - Total Supply")
print(f"Deficit = {q_total_need:.2f} - {q_total_supply:.2f} = {deficit:.2f} MVAr\n")

# Final equation summary
print("--- Final Equation with All Numbers ---")
print("Total Need = ( (Q_load1 + Q_load2 + Q_load3) * (1 + %increase) ) + Q_loss")
print(f"Total Need = ( ({Q_load1_initial} + {Q_load2_initial} + {Q_load3_initial}) * (1 + {q_demand_increase_factor}) ) + {q_loss_line:.2f} MVAr = {q_total_need:.2f} MVAr")
print("Total Supply = (Q_load1 * %S3_support) + Q_battery")
print(f"Total Supply = ({Q_load1_initial} * {s3_support_factor}) + {battery_support} = {q_total_supply:.2f} MVAr")
print("\nDeficit = Total Need - Total Supply")
print(f"Deficit = {q_total_need:.2f} MVAr - {q_total_supply:.2f} MVAr = {deficit:.2f} MVAr")