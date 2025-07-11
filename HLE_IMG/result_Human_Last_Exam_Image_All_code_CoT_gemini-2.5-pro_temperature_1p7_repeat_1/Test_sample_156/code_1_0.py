import math

# Step 1: Define initial system parameters from the problem description
# Load data (MW, MVAr)
p_load1 = 125.0
q_load1_initial = 50.0

p_load2 = 90.0
q_load2_initial = 30.0

p_load3 = 100.0
q_load3_initial = 35.0

# Transmission line data
line_reactance_per_km = 0.12  # in Ohm/km
line_length = 80.0  # in km
line_voltage = 230.0  # in kV

# Reactive power support data
s3_support_fraction = 0.50
battery_support_q = 10.0  # in MVAr

# System condition changes
q_demand_increase_factor = 1.10  # 10% increase

print("Step-by-step calculation for the reactive power deficit:\n")

# Step 2: Calculate the total updated reactive power demand
q_load1_updated = q_load1_initial * q_demand_increase_factor
q_load2_updated = q_load2_initial * q_demand_increase_factor
q_load3_updated = q_load3_initial * q_demand_increase_factor
q_total_updated = q_load1_updated + q_load2_updated + q_load3_updated

print(f"1. Updated Total Reactive Power Demand (with 10% increase):")
print(f"   Q_demand_updated = (Q_Load1_initial + Q_Load2_initial + Q_Load3_initial) * {q_demand_increase_factor}")
print(f"   Q_demand_updated = ({q_load1_initial} + {q_load2_initial} + {q_load3_initial}) * {q_demand_increase_factor} = {q_total_updated:.2f} MVAr\n")


# Step 3: Calculate reactive power losses in the transmission line
# First, find the total power flowing to calculate the loss. We assume the total load is used for this calculation.
p_total = p_load1 + p_load2 + p_load3
# Apparent power S = sqrt(P^2 + Q^2). Use the updated Q demand.
s_total_updated = math.sqrt(p_total**2 + q_total_updated**2)
# Calculate the total line reactance
line_reactance_total = line_reactance_per_km * line_length
# Calculate reactive power loss using the formula Q_loss = (|S|^2 / |V|^2) * X
# S is in MVA, V in kV, X in Ohms. The result will be in MVAr.
q_loss = (s_total_updated**2 / line_voltage**2) * line_reactance_total
print(f"2. Reactive Power Loss in the Transmission Line (Bus 7 to Bus 9):")
print(f"   Total Real Power P_total = {p_load1} + {p_load2} + {p_load3} = {p_total:.2f} MW")
print(f"   Total Updated Apparent Power S_updated = sqrt(P_total^2 + Q_demand_updated^2)")
print(f"   S_updated = sqrt({p_total:.2f}^2 + {q_total_updated:.2f}^2) = {s_total_updated:.2f} MVA")
print(f"   Total Line Reactance X_line = {line_reactance_per_km} Ohm/km * {line_length} km = {line_reactance_total:.2f} Ohms")
print(f"   Q_loss = (S_updated^2 / V_line^2) * X_line")
print(f"   Q_loss = ({s_total_updated:.2f}^2 / {line_voltage:.1f}^2) * {line_reactance_total:.2f} = {q_loss:.2f} MVAr\n")
print(f"   Note: The skin effect increases resistance by 5%, which affects real power (P) loss. It is not needed for the reactive power (Q) loss calculation.\n")

# Step 4: Calculate the total available reactive power support
# S3 provides 50% of Load 1's *initial* reactive power
q_support_s3 = q_load1_initial * s3_support_fraction
q_support_total = q_support_s3 + battery_support_q
print(f"3. Total Available Reactive Power Support:")
print(f"   Q_support_S3 = {s3_support_fraction} * Q_Load1_initial = {s3_support_fraction} * {q_load1_initial} = {q_support_s3:.2f} MVAr")
print(f"   Q_support_Battery = {battery_support_q:.2f} MVAr")
print(f"   Total Q_support = {q_support_s3:.2f} + {battery_support_q:.2f} = {q_support_total:.2f} MVAr\n")

# Step 5: Calculate the final reactive power deficit
# Deficit = (Total Demand + Losses) - Total Support
total_q_needed = q_total_updated + q_loss
deficit = total_q_needed - q_support_total
print(f"4. Final Reactive Power Deficit Calculation:")
print(f"   Total Reactive Power Needed = Q_demand_updated + Q_loss")
print(f"   Total Reactive Power Needed = {q_total_updated:.2f} + {q_loss:.2f} = {total_q_needed:.2f} MVAr")
print(f"   Final Deficit = Total Reactive Power Needed - Total Q_support")
print(f"   Final Deficit = {total_q_needed:.2f} - {q_support_total:.2f} = {deficit:.2f} MVAr")
<<<112.41>>>