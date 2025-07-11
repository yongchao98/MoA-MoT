import math

# --- Given Data ---

# Load Information (from diagram)
P_load1 = 125.0  # MW
Q_load1_initial = 50.0  # MVAr

P_load2 = 90.0  # MW
Q_load2_initial = 30.0  # MVAr

P_load3 = 100.0  # MW
Q_load3_initial = 35.0  # MVAr

# Line Impedance and Voltage
Z_per_km_imag = 0.12  # Ohm/km (Reactance)
line_length = 80.0  # km
line_voltage = 230.0  # kV

# Reactive Power Support
S3_support_factor = 0.50  # 50%
BESS_support = 10.0  # MVAr

# System Condition Factor
reactive_demand_increase = 0.10  # 10%

# --- Calculation ---

# Step 1: Calculate the total initial and adjusted reactive power demand.
Q_total_initial = Q_load1_initial + Q_load2_initial + Q_load3_initial
Q_total_adjusted = Q_total_initial * (1 + reactive_demand_increase)

# Step 2: Calculate the reactive power loss in the transmission line.
# 2a: Calculate total real power and total apparent power of the loads.
P_total = P_load1 + P_load2 + P_load3
S_total_MVA = math.sqrt(P_total**2 + Q_total_adjusted**2)

# 2b: Calculate total line reactance (X_line).
X_line = Z_per_km_imag * line_length

# 2c: Calculate reactive power loss (Q_loss) using Q_loss = X * |S|^2 / |V|^2
Q_loss = X_line * (S_total_MVA**2) / (line_voltage**2)

# Step 3: Calculate the total reactive power requirement.
Q_required = Q_total_adjusted + Q_loss

# Step 4: Calculate the total available reactive power supply.
Q_S3_support = S3_support_factor * Q_load1_initial
Q_total_supply = Q_S3_support + BESS_support

# Step 5: Calculate the final reactive power deficit.
deficit = Q_required - Q_total_supply

# --- Output the results ---
print("Calculation of Reactive Power Deficit")
print("="*40)
print(f"1. Total Adjusted Reactive Load Demand:")
print(f"   Initial Demand = {Q_load1_initial:.1f} + {Q_load2_initial:.1f} + {Q_load3_initial:.1f} = {Q_total_initial:.1f} MVAr")
print(f"   Adjusted Demand (with 10% increase) = {Q_total_initial:.1f} * (1 + {reactive_demand_increase}) = {Q_total_adjusted:.1f} MVAr\n")

print(f"2. Reactive Power Line Loss:")
print(f"   Total Apparent Power |S| = sqrt({P_total:.1f}² + {Q_total_adjusted:.1f}²) = {S_total_MVA:.2f} MVA")
print(f"   Total Line Reactance X = {Z_per_km_imag} * {line_length} = {X_line:.1f} Ω")
print(f"   Line Loss = {X_line:.1f} * {S_total_MVA:.2f}² / {line_voltage}² = {Q_loss:.2f} MVAr\n")

print(f"3. Total Reactive Power Requirement:")
print(f"   Requirement = Adjusted Demand + Line Loss")
print(f"   Requirement = {Q_total_adjusted:.1f} + {Q_loss:.2f} = {Q_required:.2f} MVAr\n")

print(f"4. Total Available Reactive Power Supply:")
print(f"   Supply from S3 = {S3_support_factor} * {Q_load1_initial} = {Q_S3_support:.1f} MVAr")
print(f"   Supply from BESS = {BESS_support:.1f} MVAr")
print(f"   Total Supply = {Q_S3_support:.1f} + {BESS_support:.1f} = {Q_total_supply:.1f} MVAr\n")

print(f"5. Final Reactive Power Deficit:")
print(f"   Deficit = Total Requirement - Total Supply")
print(f"   Deficit = {Q_required:.2f} MVAr - {Q_total_supply:.1f} MVAr")
print(f"   Final Deficit = {deficit:.2f} MVAr\n")

# Final equation with all numbers
print("Final Equation:")
print(f"Deficit = (({Q_load1_initial} + {Q_load2_initial} + {Q_load3_initial}) * (1 + {reactive_demand_increase}) + {Q_loss:.2f}) - ({S3_support_factor} * {Q_load1_initial} + {BESS_support})")
print(f"Deficit = ({Q_total_adjusted:.1f} + {Q_loss:.2f}) - ({Q_S3_support:.1f} + {BESS_support:.1f})")
print(f"Deficit = {Q_required:.2f} - {Q_total_supply:.1f} = {deficit:.2f} MVAr")