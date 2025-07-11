import math

# --- Step 1: Calculate Total Adjusted Reactive Power Demand ---

# Initial reactive power demands from the diagram (in MVAr)
q_load1_initial = 50.0
q_load2_initial = 30.0
q_load3_initial = 35.0

# Increase factor due to voltage drop
demand_increase_factor = 1.10

# Calculate adjusted reactive power demand for each load
q_load1_adj = q_load1_initial * demand_increase_factor
q_load2_adj = q_load2_initial * demand_increase_factor
q_load3_adj = q_load3_initial * demand_increase_factor

# Calculate total adjusted reactive power demand
total_q_demand = q_load1_adj + q_load2_adj + q_load3_adj

print(f"Step 1: Calculating Total Adjusted Reactive Power Demand")
print(f"Initial Reactive Power for Load 1: {q_load1_initial} MVAr")
print(f"Initial Reactive Power for Load 2: {q_load2_initial} MVAr")
print(f"Initial Reactive Power for Load 3: {q_load3_initial} MVAr")
print(f"After a 10% increase, the total reactive power demand is: {total_q_demand:.2f} MVAr\n")


# --- Step 2: Calculate Reactive Power Loss ---

# Line parameters
line_impedance_reactive_per_km = 0.12  # in Ohm/km
line_length_km = 80.0
line_voltage_kv = 230.0

# Real power demands from the diagram (in MW)
p_load1 = 125.0
p_load2 = 90.0
p_load3 = 100.0
total_p_demand = p_load1 + p_load2 + p_load3

# Assumption: Calculate losses based on total system load flowing through the line.
# Calculate total apparent power (S) of the system
total_s_demand_mva = math.sqrt(total_p_demand**2 + total_q_demand**2)

# Calculate total line reactance (X)
line_reactance_ohm = line_impedance_reactive_per_km * line_length_km

# Calculate reactive power loss (Q_loss) using the formula: Q_loss = (|S|^2 / |V|^2) * X
# S must be in MVA, V in kV for the result to be in MVAr
q_loss_mvar = (total_s_demand_mva**2 / line_voltage_kv**2) * line_reactance_ohm

print(f"Step 2: Calculating Reactive Power Loss in the Transmission Line")
print(f"Total Real Power Demand: {total_p_demand} MW")
print(f"Total Apparent Power Demand (S): {total_s_demand_mva:.2f} MVA")
print(f"Total Line Reactance (X): {line_reactance_ohm:.2f} Ohms")
print(f"Reactive Power Loss (Q_loss): {q_loss_mvar:.2f} MVAr\n")


# --- Step 3: Calculate Total Available Reactive Power ---

# Reactive power supplied by Generator S3 (50% of Load 1's initial Q)
q_supply_s3 = 0.50 * q_load1_initial

# Reactive power supplied by the Battery (BESS)
q_supply_bess = 10.0

# Total available reactive power supply
total_q_supply = q_supply_s3 + q_supply_bess

print(f"Step 3: Calculating Available Reactive Power Supply")
print(f"Reactive Power from Generator S3: {q_supply_s3:.2f} MVAr")
print(f"Reactive Power from BESS: {q_supply_bess:.2f} MVAr")
print(f"Total Available Supply: {total_q_supply:.2f} MVAr\n")


# --- Step 4: Calculate the Reactive Power Deficit ---

# Total reactive power required is the sum of total demand and losses
total_q_required = total_q_demand + q_loss_mvar

# Deficit is the difference between required and supplied reactive power
q_deficit = total_q_required - total_q_supply

print(f"Step 4: Calculating the Final Reactive Power Deficit")
print(f"Total Reactive Power Required (Demand + Loss): {total_q_required:.2f} MVAr")
print(f"Total Reactive Power Supplied: {total_q_supply:.2f} MVAr")
print(f"Reactive Power Deficit = Required - Supplied = {q_deficit:.2f} MVAr\n")

# --- Final Equation Output ---
print("Final Calculation:")
print(f"Deficit = (Q_Load1_adj + Q_Load2_adj + Q_Load3_adj + Q_Loss) - (Q_S3_Supply + Q_BESS_Supply)")
print(f"Deficit = ({q_load1_adj:.2f} + {q_load2_adj:.2f} + {q_load3_adj:.2f} + {q_loss_mvar:.2f}) - ({q_supply_s3:.2f} + {q_supply_bess:.2f}) = {q_deficit:.2f} MVAr")
