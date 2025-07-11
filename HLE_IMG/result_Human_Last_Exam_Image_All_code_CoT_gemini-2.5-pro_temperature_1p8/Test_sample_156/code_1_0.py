import math

# Step 0: Define initial parameters from the problem
# Load Information (Initial)
p_load1 = 125  # MW
q_load1_initial = 50  # MVAr

p_load2 = 90   # MW
q_load2_initial = 30  # MVAr

p_load3 = 100  # MW
q_load3_initial = 35  # MVAr

# Line Information
z_real_per_km = 0.05  # Ohm/km
z_imag_per_km = 0.12  # Ohm/km
line_length = 80      # km
line_voltage = 230    # kV

# System Adjustments and Support
q_increase_factor = 1.10  # 10% increase
s3_support_factor = 0.50 # 50%
q_bess_support = 10       # MVAr

# --- Calculations ---

# Step 1: Calculate the final reactive power demand from loads
print("Step 1: Calculate the final reactive power demand from loads after the 10% increase.")
q_load1_final = q_load1_initial * q_increase_factor
q_load2_final = q_load2_initial * q_increase_factor
q_load3_final = q_load3_initial * q_increase_factor
total_q_loads = q_load1_final + q_load2_final + q_load3_final
print(f"Final reactive demand for Load 1: {q_load1_initial} MVAr * {q_increase_factor} = {q_load1_final:.2f} MVAr")
print(f"Final reactive demand for Load 2: {q_load2_initial} MVAr * {q_increase_factor} = {q_load2_final:.2f} MVAr")
print(f"Final reactive demand for Load 3: {q_load3_initial} MVAr * {q_increase_factor} = {q_load3_final:.2f} MVAr")
print(f"Total load reactive power demand = {q_load1_final:.2f} + {q_load2_final:.2f} + {q_load3_final:.2f} = {total_q_loads:.2f} MVAr\n")


# Step 2: Calculate reactive power losses in the transmission line
print("Step 2: Calculate reactive power losses in the transmission line.")
# Calculate total power to find the line current (assuming total load flows through the line)
total_p_loads = p_load1 + p_load2 + p_load3
total_s_loads_mva = math.sqrt(total_p_loads**2 + total_q_loads**2)
print(f"Total System Real Power (P): {total_p_loads:.2f} MW")
print(f"Total System Reactive Power (Q): {total_q_loads:.2f} MVAr")
print(f"Total System Apparent Power (S): sqrt({total_p_loads:.2f}^2 + {total_q_loads:.2f}^2) = {total_s_loads_mva:.2f} MVA")

# Calculate line current
line_current_a = (total_s_loads_mva * 1e6) / (math.sqrt(3) * line_voltage * 1e3)
print(f"Line Current (I): {total_s_loads_mva:.2f} MVA / (sqrt(3) * {line_voltage} kV) = {line_current_a:.2f} A")

# Calculate total line reactance
x_total = z_imag_per_km * line_length
print(f"Total Line Reactance (X): {z_imag_per_km} Ohm/km * {line_length} km = {x_total:.2f} Ohms")

# Calculate reactive power loss (Q = 3 * I^2 * X for three-phase)
q_loss_var = 3 * (line_current_a**2) * x_total
q_loss_mvar = q_loss_var / 1e6
print(f"Total line reactive power loss = 3 * ({line_current_a:.2f} A)^2 * {x_total:.2f} Ohms = {q_loss_mvar:.2f} MVAr\n")


# Step 3: Calculate total reactive power demand (Loads + Losses)
print("Step 3: Calculate total reactive power demand.")
total_q_demand = total_q_loads + q_loss_mvar
print(f"Total Reactive Power Demand = Load Demand + Line Losses")
print(f"Total Reactive Power Demand = {total_q_loads:.2f} MVAr + {q_loss_mvar:.2f} MVAr = {total_q_demand:.2f} MVAr\n")


# Step 4: Calculate total available reactive power supply
print("Step 4: Calculate total available reactive power supply.")
q_s3_support = q_load1_initial * s3_support_factor
total_q_supply = q_s3_support + q_bess_support
print(f"Reactive Power from S3: {s3_support_factor:.2f} * {q_load1_initial} MVAr (Load 1 initial) = {q_s3_support:.2f} MVAr")
print(f"Reactive Power from BESS: {q_bess_support:.2f} MVAr")
print(f"Total Reactive Power Supply = {q_s3_support:.2f} MVAr + {q_bess_support:.2f} MVAr = {total_q_supply:.2f} MVAr\n")


# Step 5: Calculate the total reactive power deficit
print("Step 5: Calculate the final total reactive power deficit.")
q_deficit = total_q_demand - total_q_supply
print(f"Total Reactive Power Deficit = Total Demand - Total Supply")
print(f"Total Reactive Power Deficit = {total_q_demand:.2f} MVAr - {total_q_supply:.2f} MVAr = {q_deficit:.2f} MVAr\n")

print("The equation for the final answer is:")
print(f"Deficit = (Q_Load1_final + Q_Load2_final + Q_Load3_final + Q_Line_Loss) - (Q_S3_Support + Q_BESS_Support)")
print(f"Deficit = ({q_load1_final:.2f} + {q_load2_final:.2f} + {q_load3_final:.2f} + {q_loss_mvar:.2f}) - ({q_s3_support:.2f} + {q_bess_support:.2f}) = {q_deficit:.2f} MVAr")

# Final Answer
# print(f"<<<{q_deficit:.2f}>>>")