import cmath

# --- Step 1: Initial Reactive Power Demand ---
# From the diagram
q_load1_initial = 50  # MVAr
q_load2_initial = 30  # MVAr
q_load3_initial = 35  # MVAr

# Calculate initial total demand
total_q_demand_initial = q_load1_initial + q_load2_initial + q_load3_initial

print("Step 1: Calculate Initial Total Reactive Power Demand")
print(f"Reactive Power at Load 1: {q_load1_initial} MVAr")
print(f"Reactive Power at Load 2: {q_load2_initial} MVAr")
print(f"Reactive Power at Load 3: {q_load3_initial} MVAr")
print(f"Initial Total Reactive Demand = {q_load1_initial} + {q_load2_initial} + {q_load3_initial} = {total_q_demand_initial} MVAr\n")

# --- Step 2: Increased Reactive Power Demand ---
# A 10% increase occurs at each load due to voltage drop
increase_percentage = 10
increase_factor = 1 + (increase_percentage / 100)

# Calculate final reactive power for each load
q_load1_final = q_load1_initial * increase_factor
q_load2_final = q_load2_initial * increase_factor
q_load3_final = q_load3_initial * increase_factor

# Calculate final total demand
total_q_demand_final = q_load1_final + q_load2_final + q_load3_final

print(f"Step 2: Calculate Final Total Reactive Power Demand with {increase_percentage}% Increase")
print(f"New Reactive Power at Load 1 = {q_load1_initial} * {increase_factor} = {q_load1_final} MVAr")
print(f"New Reactive Power at Load 2 = {q_load2_initial} * {increase_factor} = {q_load2_final} MVAr")
print(f"New Reactive Power at Load 3 = {q_load3_initial} * {increase_factor:.1f} = {q_load3_final:.1f} MVAr")
print(f"Final Total Reactive Demand = {q_load1_final} + {q_load2_final} + {q_load3_final:.1f} = {total_q_demand_final:.1f} MVAr\n")

# --- Step 3: Available Reactive Power Supply ---
# Generator S3 provides 50% of Load 1's reactive power
s3_support_percentage = 50
q_s3 = (s3_support_percentage / 100) * q_load1_initial

# BESS provides a fixed amount
q_bess = 10  # MVAr

# Calculate total supply
total_q_supply = q_s3 + q_bess

print("Step 3: Calculate Total Available Reactive Power Supply")
print(f"Reactive Power from Generator S3 ({s3_support_percentage}% of Load 1) = {s3_support_percentage/100} * {q_load1_initial} = {q_s3} MVAr")
print(f"Reactive Power from BESS = {q_bess} MVAr")
print(f"Total Reactive Supply = {q_s3} + {q_bess} = {total_q_supply} MVAr\n")

# --- Step 4: Reactive Power Deficit ---
# Deficit = Demand - Supply
deficit = total_q_demand_final - total_q_supply

print("Step 4: Calculate the Reactive Power Deficit")
print("Reactive Power Deficit = Final Total Demand - Total Supply")
print(f"Reactive Power Deficit = {total_q_demand_final:.1f} - {total_q_supply} = {deficit:.1f} MVAr")

# Final Answer
# print(f"\nFinal Answer: {deficit:.1f}")
<<<91.5>>>