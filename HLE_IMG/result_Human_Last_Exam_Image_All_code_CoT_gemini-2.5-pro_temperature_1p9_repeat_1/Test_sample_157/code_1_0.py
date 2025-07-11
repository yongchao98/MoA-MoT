import math

# Step 1: Define system parameters from the problem description
# Load data
base_load_per_bus = 50.0 / 3.0  # MW
num_base_loads = 3
additional_load = 100.0  # MW, problem specifies the *active* load is connected

# System operational parameters
# We interpret "system capacity" as the transformer rating of 210 MVA.
# This is a reasonable assumption as using the sum of generator ratings (300 MVA)
# results in an excessively high power loss.
system_capacity_S = 210.0  # MVA
power_factor = 0.9

# Step 2: Calculate the total active power load (P_load)
total_base_load = num_base_loads * base_load_per_bus
total_load_P = total_base_load + additional_load

# Step 3: Calculate the total generated active power (P_gen)
total_gen_P = system_capacity_S * power_factor

# Step 4: Calculate the total power loss (P_loss)
# P_loss = P_gen - P_load
total_loss_P = total_gen_P - total_load_P

# Step 5: Print the results following the calculation steps.
print("Step 1: Calculate Total Active Power Load (Pl)")
print(f"Pl = (Number of Base Loads * Power per Base Load) + Additional Load")
print(f"Pl = ({num_base_loads} * {base_load_per_bus:.3f} MW) + {additional_load:.3f} MW = {total_load_P:.3f} MW\n")

print("Step 2: Calculate Total Generated Active Power (Pg)")
print(f"Pg = System Capacity * Power Factor")
print(f"Pg = {system_capacity_S:.3f} MVA * {power_factor:.3f} = {total_gen_P:.3f} MW\n")

print("Step 3: Calculate Total Power Loss")
print(f"Total Power Loss = Pg - Pl")
print(f"Total Power Loss = {total_gen_P:.3f} MW - {total_load_P:.3f} MW = {total_loss_P:.3f} MW")
