import math

# --- Step 1: Calculate Total Active Power Generation ---

# Apparent power rating of each generator (MVA)
S_gen_unit = 100.0

# Number of generators
num_generators = 3

# Power factor
power_factor = 0.9

# Total apparent power generation (MVA)
S_gen_total = S_gen_unit * num_generators

# Total active power generation (MW)
P_gen_total = S_gen_total * power_factor

# --- Step 2: Calculate Total Active Power Load ---

# Base load at each of the three buses (MW)
P_base_load_unit = 50.0 / 3.0

# Number of base loads
num_base_loads = 3

# Total base load (MW)
P_base_load_total = P_base_load_unit * num_base_loads

# Additional active load at Bus 8 (MW)
P_additional_load = 100.0

# Total active power load (MW)
P_load_total = P_base_load_total + P_additional_load

# --- Step 3: Calculate Total Power Loss ---

# Total power loss is the difference between generation and load (MW)
P_loss_total = P_gen_total - P_load_total

# --- Step 4: Print the results ---

print("This script calculates the total resistive power losses in the system.")
print("\nStep 1: Calculate Total Active Power Generation")
print(f"Total Apparent Power Generation (S_gen) = {num_generators} generators * {S_gen_unit:.3f} MVA/generator = {S_gen_total:.3f} MVA")
print(f"Total Active Power Generation (P_gen) = S_gen * Power Factor = {S_gen_total:.3f} MVA * {power_factor:.3f} = {P_gen_total:.3f} MW")

print("\nStep 2: Calculate Total Active Power Load")
print(f"Total Base Load = {num_base_loads} loads * {P_base_load_unit:.3f} MW/load = {P_base_load_total:.3f} MW")
print(f"Additional Load = {P_additional_load:.3f} MW")
print(f"Total Active Load (P_load) = Total Base Load + Additional Load = {P_base_load_total:.3f} MW + {P_additional_load:.3f} MW = {P_load_total:.3f} MW")

print("\nStep 3: Calculate Total Power Loss")
print("The total power loss is the difference between total active power generation and total active load.")
print("Final Equation:")
print(f"Total Power Loss (MW) = P_gen (MW) - P_load (MW)")
print(f"Total Power Loss (MW) = {P_gen_total:.3f} - {P_load_total:.3f}")
print(f"Total Power Loss (MW) = {P_loss_total:.3f}")

# The final answer in the required format
# print(f"\n<<<{P_loss_total:.3f}>>>")