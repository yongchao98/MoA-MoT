import math

# Step 1: Define initial parameters from the problem description
# Load Information (initial)
p_load1 = 125  # MW
q_load1_initial = 50  # MVAr

p_load2 = 90  # MW
q_load2_initial = 30  # MVAr

p_load3 = 100  # MW
q_load3_initial = 35  # MVAr

# Line Information
z_line_reactance_per_km = 0.12  # Ω/km
line_length = 80  # km
line_voltage = 230  # kV

# Voltage Drop Effect
demand_increase_factor = 1.10  # 10% increase

# Reactive Power Support
s3_support_factor = 0.50  # 50% of Load 1's reactive power
bess_support = 10  # MVAr

# Step 2: Calculate the final reactive power demand of the loads
q_load1_final = q_load1_initial * demand_increase_factor
q_load2_final = q_load2_initial * demand_increase_factor
q_load3_final = q_load3_initial * demand_increase_factor
q_demand_total = q_load1_final + q_load2_final + q_load3_final

# Step 3: Calculate the total real power demand
p_demand_total = p_load1 + p_load2 + p_load3

# Step 4: Estimate the "full load" current to calculate line losses
# Calculate total apparent power (S) from total real (P) and reactive (Q) power
s_demand_total = math.sqrt(p_demand_total**2 + q_demand_total**2)

# Calculate current in kA: I = S(MVA) / (sqrt(3) * V(kV))
# Note: The skin effect on resistance is not needed for reactive power loss calculation.
current_kA = s_demand_total / (math.sqrt(3) * line_voltage)

# Calculate total line reactance
x_line_total = z_line_reactance_per_km * line_length

# Calculate reactive power loss in MVAr: Q_loss = I(kA)^2 * X(Ω)
q_loss = (current_kA**2) * x_line_total

# Step 5: Calculate the total reactive power requirement
q_requirement = q_demand_total + q_loss

# Step 6: Calculate the total specified reactive power supply
# S3 support is based on the initial reactive power of Load 1
q_s3_support = s3_support_factor * q_load1_initial
q_supply_total = q_s3_support + bess_support

# Step 7: Calculate the reactive power deficit
deficit = q_requirement - q_supply_total

# Print the results step-by-step
print("--- Calculation of Reactive Power Deficit ---")
print("\n1. Total Reactive Power Demand (after 10% increase):")
print(f"   Increased Demand = ({q_load1_initial} + {q_load2_initial} + {q_load3_initial}) * {demand_increase_factor} = {q_demand_total:.2f} MVAr")

print("\n2. Reactive Power Loss in Transmission Line (Bus 7 to 9):")
print(f"   Total Line Reactance = {z_line_reactance_per_km} Ω/km * {line_length} km = {x_line_total:.2f} Ω")
print(f"   Estimated Full Load Current = {current_kA:.2f} kA")
print(f"   Reactive Power Loss = ({current_kA:.2f})^2 * {x_line_total:.2f} = {q_loss:.2f} MVAr")

print("\n3. Total Reactive Power Requirement:")
print(f"   Total Requirement = Demand + Loss = {q_demand_total:.2f} + {q_loss:.2f} = {q_requirement:.2f} MVAr")

print("\n4. Available Reactive Power Support:")
print(f"   S3 Support = {s3_support_factor} * {q_load1_initial} = {q_s3_support:.2f} MVAr")
print(f"   BESS Support = {bess_support:.2f} MVAr")
print(f"   Total Support = {q_s3_support:.2f} + {bess_support:.2f} = {q_supply_total:.2f} MVAr")

print("\n5. Final Reactive Power Deficit:")
print(f"   Deficit = Total Requirement - Total Support")
print(f"   Deficit = {q_requirement:.2f} MVAr - {q_supply_total:.2f} MVAr = {deficit:.2f} MVAr")

print("\nThe total reactive power deficit is {:.2f} MVAr.".format(deficit))
print("<<<98.48>>>")