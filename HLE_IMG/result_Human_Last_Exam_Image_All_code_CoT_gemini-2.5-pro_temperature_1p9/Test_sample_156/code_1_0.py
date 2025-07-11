import cmath

# 1. Define initial parameters from the problem description
# Load Information (initial values)
P_load1 = 125  # MW
Q_load1_initial = 50  # MVAr

P_load2 = 90   # MW
Q_load2_initial = 30   # MVAr

P_load3 = 100  # MW
Q_load3_initial = 35  # MVAr

# Voltage Drop Effect
q_increase_factor = 1.10

# Reactive Power Support
Q_S3_support_factor = 0.50
Q_BESS_support = 10  # MVAr

# Line Information
line_length_total = 80  # km
line_reactance_per_km = 0.12  # Ω/km
line_voltage = 230e3  # V

# --- Step-by-step Calculation ---

# Step 1: Calculate adjusted reactive power demand from loads
Q_load1_adj = Q_load1_initial * q_increase_factor
Q_load2_adj = Q_load2_initial * q_increase_factor
Q_load3_adj = Q_load3_initial * q_increase_factor
Q_loads_total_adj = Q_load1_adj + Q_load2_adj + Q_load3_adj

# Step 2: Calculate reactive power supply from specified sources
Q_S3_supply = Q_load1_initial * Q_S3_support_factor
Q_supply_specified = Q_S3_supply + Q_BESS_support

# Step 3: Calculate reactive power loss in the transmission line
# Assumption: Bus 8 is at the midpoint of the 80km line.
length_78 = line_length_total / 2
length_89 = line_length_total / 2

# Calculate reactance for each line segment
X_78 = line_reactance_per_km * length_78
X_89 = line_reactance_per_km * length_89

# Calculate complex power flows (S = P + jQ) in MVA
# Net power drawn from network at Bus 9 (Load 2 demand - local generation)
S_L2_adj = complex(P_load2, Q_load2_adj)
Q_gen_at_B9 = Q_S3_supply + Q_BESS_support
S_net_draw_at_B9 = S_L2_adj - complex(0, Q_gen_at_B9)
S_flow_89 = S_net_draw_at_B9

# Net power drawn at Bus 8
S_L3_adj = complex(P_load3, Q_load3_adj)

# Power flow from Bus 7 to Bus 8
S_flow_78 = S_flow_89 + S_L3_adj

# Calculate reactive power loss (Q_loss = |S|^2 * X / |V|^2)
# Convert MVA to VA for calculation
S_flow_78_va = S_flow_78 * 1e6
S_flow_89_va = S_flow_89 * 1e6

Q_loss_78 = (abs(S_flow_78_va)**2 * X_78) / (line_voltage**2) / 1e6 # in MVAr
Q_loss_89 = (abs(S_flow_89_va)**2 * X_89) / (line_voltage**2) / 1e6 # in MVAr

Q_loss_line_total = Q_loss_78 + Q_loss_89

# Step 4: Calculate total demand and the final deficit
Q_demand_total = Q_loads_total_adj + Q_loss_line_total
Q_deficit = Q_demand_total - Q_supply_specified


# --- Output the results ---
print("--- Calculation of Reactive Power Deficit ---")
print("\nThe final deficit is calculated as: Total Demand - Specified Supply.")
print("Total Demand = (Adjusted Load Demands) + (Line Losses)")
print("Specified Supply = (Generator S3 Support) + (BESS Support)")
print("\nHere are the components of the equation:\n")

print(f"Adjusted Reactive Power Demand for Load 1: {Q_load1_adj:.2f} MVAr")
print(f"Adjusted Reactive Power Demand for Load 2: {Q_load2_adj:.2f} MVAr")
print(f"Adjusted Reactive Power Demand for Load 3: {Q_load3_adj:.2f} MVAr")
print(f"Total Reactive Power Loss in Transmission Line: {Q_loss_line_total:.2f} MVAr")
print(f"Reactive Power Supplied by Generator S3: {Q_S3_supply:.2f} MVAr")
print(f"Reactive Power Supplied by BESS: {Q_BESS_support:.2f} MVAr")

print("\n--- Final Equation ---")
print(f"Total Reactive Demand = {Q_load1_adj:.2f} + {Q_load2_adj:.2f} + {Q_load3_adj:.2f} + {Q_loss_line_total:.2f} = {Q_demand_total:.2f} MVAr")
print(f"Total Specified Supply = {Q_S3_supply:.2f} + {Q_BESS_support:.2f} = {Q_supply_specified:.2f} MVAr")
print(f"Reactive Power Deficit = {Q_demand_total:.2f} - {Q_supply_specified:.2f} = {Q_deficit:.2f} MVAr")

print("\nThe total reactive power deficit is:")
print(f"{Q_deficit:.2f} MVAr")

# Final answer in specified format
print(f"\n<<<{Q_deficit:.2f}>>>")