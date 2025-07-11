import math

# --- Given Data ---

# Load Information (from diagram)
P_load1 = 125.0  # MW
Q_load1_initial = 50.0  # MVAr

P_load2 = 90.0  # MW
Q_load2_initial = 30.0  # MVAr

P_load3 = 100.0  # MW
Q_load3_initial = 35.0  # MVAr

# Line Impedance and characteristics
z_per_km_real = 0.05  # Ohm/km
z_per_km_imag = 0.12  # Ohm/km
line_length = 80.0  # km
line_voltage = 230.0  # kV

# System Factors
skin_effect_increase = 0.05  # 5%
q_demand_increase = 0.10  # 10%

# Reactive Power Support
q_s3_support_factor = 0.50  # 50% of Load 1's reactive power
q_bess_support = 10.0  # MVAr

# --- Calculations ---

# Step 1: Calculate the new reactive power demand for each load
q_load1_new = Q_load1_initial * (1 + q_demand_increase)
q_load2_new = Q_load2_initial * (1 + q_demand_increase)
q_load3_new = Q_load3_initial * (1 + q_demand_increase)

total_q_demand_at_loads = q_load1_new + q_load2_new + q_load3_new

# Step 2: Calculate reactive power loss in the transmission line (Bus 7 to Bus 9)
# Assumption: Bus 8 is at the midpoint of the line.
# The line is considered in two segments: Bus 7 to Bus 8 and Bus 8 to Bus 9.

# Impedance for one segment (40 km)
segment_length = line_length / 2.0
x_segment = z_per_km_imag * segment_length

# Complex power (S = P + jQ) for loads with increased reactive demand
S_load2_new = complex(P_load2, q_load2_new)
S_load3_new = complex(P_load3, q_load3_new)

# Reactive power support at the Bus 9 end
q_s3_support = Q_load1_initial * q_s3_support_factor
total_q_support_at_9 = q_s3_support + q_bess_support
S_support_at_9 = complex(0, total_q_support_at_9)

# Power flow from Bus 8 to Bus 9
# This flow must supply Load 2, minus the support available at Bus 9
S_flow_8_9 = S_load2_new - S_support_at_9
S_flow_8_9_mag_MVA = abs(S_flow_8_9)

# Reactive power loss in segment 8-9 using Q_loss = X * |S|^2 / |V|^2
# Result is converted from VAr to MVAr
q_loss_8_9 = (x_segment * (S_flow_8_9_mag_MVA * 1e6)**2) / ((line_voltage * 1e3)**2) / 1e6

# Power flow from Bus 7 to Bus 8
# This flow must supply Load 3 plus the power flowing to the 8-9 segment
S_flow_7_8 = S_load3_new + S_flow_8_9
S_flow_7_8_mag_MVA = abs(S_flow_7_8)

# Reactive power loss in segment 7-8
q_loss_7_8 = (x_segment * (S_flow_7_8_mag_MVA * 1e6)**2) / ((line_voltage * 1e3)**2) / 1e6

# Total reactive power loss in the line
total_q_loss = q_loss_7_8 + q_loss_8_9

# Step 3: Calculate total reactive power requirement
total_q_requirement = total_q_demand_at_loads + total_q_loss

# Step 4: Calculate total available dedicated reactive power support
total_q_support = q_s3_support + q_bess_support

# Step 5: Calculate the final reactive power deficit
deficit = total_q_requirement - total_q_support

# --- Output the final equation and result ---
print("The total reactive power deficit is calculated as follows:")
print("Deficit = (New Q_Load1 + New Q_Load2 + New Q_Load3 + Q_Loss_Total) - (Q_S3 + Q_BESS)")
print("\nSubstituting the calculated values into the equation:")
print(f"Deficit = ({q_load1_new:.2f} MVAr + {q_load2_new:.2f} MVAr + {q_load3_new:.2f} MVAr + {total_q_loss:.2f} MVAr) - ({q_s3_support:.2f} MVAr + {q_bess_support:.2f} MVAr)")
print(f"Deficit = ({total_q_demand_at_loads:.2f} MVAr + {total_q_loss:.2f} MVAr) - ({total_q_support:.2f} MVAr)")
print(f"Deficit = {total_q_requirement:.2f} MVAr - {total_q_support:.2f} MVAr")
print(f"\nFinal Reactive Power Deficit = {deficit:.2f} MVAr")
