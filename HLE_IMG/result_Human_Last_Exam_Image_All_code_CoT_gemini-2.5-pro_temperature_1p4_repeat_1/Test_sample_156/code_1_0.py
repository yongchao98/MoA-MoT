import math

# Step 1: Define initial power system parameters from the problem description
# Load Information (Initial Values)
P_load1 = 125  # MW
Q_load1 = 50   # MVAr

P_load2 = 90   # MW
Q_load2 = 30   # MVAr

P_load3 = 100  # MW
Q_load3 = 35   # MVAr

# Line Impedance and characteristics
Z_per_km_imag = 0.12  # Reactive component of impedance in Ohm/km
line_length = 80      # km
line_voltage = 230    # kV

# Factors and Support
Q_demand_increase_factor = 0.10     # 10% increase in reactive power demand
S3_support_factor = 0.50            # S3 provides 50% of Load 1's Q
BESS_support = 10                   # MVAr from Battery Energy Storage System

print("--- Calculation of Reactive Power Deficit ---\n")

# Step 2: Calculate new reactive power demand for each load after the 10% increase
print("1. Calculating New Load Reactive Power Demand:")
Q_load1_new = Q_load1 * (1 + Q_demand_increase_factor)
Q_load2_new = Q_load2 * (1 + Q_demand_increase_factor)
Q_load3_new = Q_load3 * (1 + Q_demand_increase_factor)
total_Q_demand_new = Q_load1_new + Q_load2_new + Q_load3_new
print(f"   - Adjusted Load 1 Demand: {Q_load1:.1f} MVAr * {1 + Q_demand_increase_factor} = {Q_load1_new:.1f} MVAr")
print(f"   - Adjusted Load 2 Demand: {Q_load2:.1f} MVAr * {1 + Q_demand_increase_factor} = {Q_load2_new:.1f} MVAr")
print(f"   - Adjusted Load 3 Demand: {Q_load3:.1f} MVAr * {1 + Q_demand_increase_factor} = {Q_load3_new:.1f} MVAr")
print(f"   - Total Adjusted Reactive Load Demand = {total_Q_demand_new:.1f} MVAr\n")

# Step 3: Calculate total available local reactive power support
print("2. Calculating Available Local Reactive Power Support:")
Q_support_S3 = Q_load1 * S3_support_factor
total_Q_support = Q_support_S3 + BESS_support
print(f"   - Support from Generator S3: {Q_load1} MVAr * {S3_support_factor} = {Q_support_S3:.1f} MVAr")
print(f"   - Support from BESS: {BESS_support:.1f} MVAr")
print(f"   - Total Local Support = {Q_support_S3:.1f} MVAr + {BESS_support:.1f} MVAr = {total_Q_support:.1f} MVAr\n")

# Step 4: Calculate reactive power loss in the transmission line
print("3. Calculating Reactive Power Loss in the Transmission Line:")
# Calculate total real power demand
total_P_demand = P_load1 + P_load2 + P_load3
# Calculate net power to be transmitted (loads minus local support)
P_transmit = total_P_demand
Q_transmit = total_Q_demand_new - total_Q_support
# Calculate magnitude squared of apparent power transmitted
S_transmit_mag_sq = P_transmit**2 + Q_transmit**2
# Calculate total line reactance
X_line = Z_per_km_imag * line_length
# Calculate reactive power loss using Q_loss = (|S|^2 * X) / |V|^2
Q_loss = (S_transmit_mag_sq * X_line) / (line_voltage**2)

print(f"   - Net Apparent Power to Transmit (S) = {P_transmit} MW + j{Q_transmit:.1f} MVAr")
print(f"   - Total Line Reactance (X) = {Z_per_km_imag} Ohm/km * {line_length} km = {X_line:.1f} Ohms")
print(f"   - Reactive Power Loss (Q_loss) = (|S|^2 * X) / |V|^2")
print(f"   - Q_loss = ({S_transmit_mag_sq:.2f} MVA^2 * {X_line:.1f} Ohm) / ({line_voltage} kV)^2 = {Q_loss:.3f} MVAr\n")

# Step 5: Calculate the final total reactive power deficit
print("4. Calculating Final Reactive Power Deficit:")
# The deficit is the total requirement (new loads + losses) minus the local support
deficit = (total_Q_demand_new + Q_loss) - total_Q_support

# Print the final calculation breakdown
print("\n--- Final Equation ---")
print(f"Deficit = (New_Q_Load1 + New_Q_Load2 + New_Q_Load3 + Q_loss) - (Q_S3 + Q_BESS)")
print(f"Deficit = ({Q_load1_new:.1f} + {Q_load2_new:.1f} + {Q_load3_new:.1f} + {Q_loss:.3f}) - ({Q_support_S3:.1f} + {BESS_support:.1f})")
print(f"Deficit = ({total_Q_demand_new:.1f} + {Q_loss:.3f}) - {total_Q_support:.1f}")
print(f"Deficit = {total_Q_demand_new + Q_loss:.3f} - {total_Q_support:.1f}")
print(f"Total Reactive Power Deficit = {deficit:.3f} MVAr")

<<<111.026>>>