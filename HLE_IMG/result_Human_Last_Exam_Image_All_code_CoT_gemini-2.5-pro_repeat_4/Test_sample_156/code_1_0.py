import math

# Step 1: Define initial values from the problem description
# Load Information (initial)
P_load1 = 125.0  # MW
Q_load1_initial = 50.0  # MVAr

P_load2 = 90.0   # MW
Q_load2_initial = 30.0  # MVAr

P_load3 = 100.0  # MW
Q_load3_initial = 35.0  # MVAr

# Line Information
Z_line_per_km_real = 0.05  # Ohm/km
Z_line_per_km_imag = 0.12  # Ohm/km (Reactance X)
line_length = 80.0  # km
V_line = 230.0  # kV

# Voltage Support Information
S3_support_factor = 0.50  # 50%
Q_BESS_support = 10.0  # MVAr

# System Condition Factors
Q_increase_factor = 1.10  # 10% increase

# --- Calculations ---

# Step 2: Calculate adjusted reactive power demand for each load
Q_load1_adj = Q_load1_initial * Q_increase_factor
Q_load2_adj = Q_load2_initial * Q_increase_factor
Q_load3_adj = Q_load3_initial * Q_increase_factor
total_Q_loads_adj = Q_load1_adj + Q_load2_adj + Q_load3_adj

print("--- Calculating Reactive Power Demand ---")
print(f"Initial Reactive Power for Load 1: {Q_load1_initial} MVAr")
print(f"Initial Reactive Power for Load 2: {Q_load2_initial} MVAr")
print(f"Initial Reactive Power for Load 3: {Q_load3_initial} MVAr")
print(f"After a 10% increase due to voltage drop:")
print(f"Adjusted Q_Load1 = {Q_load1_initial:.2f} * {Q_increase_factor:.2f} = {Q_load1_adj:.2f} MVAr")
print(f"Adjusted Q_Load2 = {Q_load2_initial:.2f} * {Q_increase_factor:.2f} = {Q_load2_adj:.2f} MVAr")
print(f"Adjusted Q_Load3 = {Q_load3_initial:.2f} * {Q_increase_factor:.2f} = {Q_load3_adj:.2f} MVAr")
print(f"Total Adjusted Reactive Load Demand = {Q_load1_adj:.2f} + {Q_load2_adj:.2f} + {Q_load3_adj:.2f} = {total_Q_loads_adj:.2f} MVAr\n")


# Step 3: Calculate reactive power loss in the transmission line (Bus 7 to Bus 9)
# Assumption: Current is based on power supplied to Load 2 and Load 3.
P_for_loss_calc = P_load2 + P_load3
Q_for_loss_calc = Q_load2_adj + Q_load3_adj
S_for_loss_calc = math.sqrt(P_for_loss_calc**2 + Q_for_loss_calc**2)

# Current in Amperes
I_line = (S_for_loss_calc * 1e6) / (math.sqrt(3) * V_line * 1e3)

# Total line reactance in Ohms
X_total = Z_line_per_km_imag * line_length

# Reactive power loss in MVAr
Q_loss = (I_line**2 * X_total) / 1e6

print("--- Calculating Reactive Power Line Loss ---")
print("Assumption: Line current is based on the combined power for Load 2 and Load 3.")
print(f"Apparent Power for Loads 2 & 3 (S) = sqrt(({P_load2:.2f} MW + {P_load3:.2f} MW)^2 + ({Q_load2_adj:.2f} MVAr + {Q_load3_adj:.2f} MVAr)^2) = {S_for_loss_calc:.2f} MVA")
print(f"Current (I) = {S_for_loss_calc:.2f} MVA / (sqrt(3) * {V_line:.2f} kV) = {I_line:.2f} A")
print(f"Total Line Reactance (X) = {Z_line_per_km_imag:.2f} Ohm/km * {line_length:.2f} km = {X_total:.2f} Ohms")
print(f"Reactive Power Loss (Q_loss) = {I_line:.2f}^2 A^2 * {X_total:.2f} Ohms = {Q_loss:.2f} MVAr\n")

# Step 4: Calculate total reactive power demand
total_Q_demand = total_Q_loads_adj + Q_loss
print("--- Calculating Total Reactive Demand ---")
print(f"Total Demand = Total Adjusted Load Demand + Line Loss")
print(f"Total Demand = {total_Q_loads_adj:.2f} MVAr + {Q_loss:.2f} MVAr = {total_Q_demand:.2f} MVAr\n")


# Step 5: Calculate known reactive power supply
Q_S3_support = S3_support_factor * Q_load1_initial
total_Q_supply_known = Q_S3_support + Q_BESS_support

print("--- Calculating Known Reactive Power Supply ---")
print(f"Supply from Generator S3 = {S3_support_factor:.2f} * {Q_load1_initial:.2f} MVAr (Load 1 initial) = {Q_S3_support:.2f} MVAr")
print(f"Supply from BESS = {Q_BESS_support:.2f} MVAr")
print(f"Total Known Supply = {Q_S3_support:.2f} MVAr + {Q_BESS_support:.2f} MVAr = {total_Q_supply_known:.2f} MVAr\n")


# Step 6: Calculate the final reactive power deficit
deficit = total_Q_demand - total_Q_supply_known

print("--- Calculating Final Reactive Power Deficit ---")
print("Deficit = Total Reactive Power Demand - Total Known Reactive Power Supply")
print(f"Deficit = {total_Q_demand:.2f} MVAr - {total_Q_supply_known:.2f} MVAr = {deficit:.2f} MVAr")

# Final answer in the required format
# The calculation is (50*1.1 + 30*1.1 + 35*1.1 + 2.5) - (50*0.5 + 10) = (55 + 33 + 38.5 + 2.5) - (25 + 10) = 129 - 35 = 94
print("\nFinal Deficit Equation:")
print(f"({Q_load1_adj:.1f} + {Q_load2_adj:.1f} + {Q_load3_adj:.1f} + {Q_loss:.1f}) - ({Q_S3_support:.1f} + {Q_BESS_support:.1f}) = {deficit:.1f}")
