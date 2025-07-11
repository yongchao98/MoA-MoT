import math

# --- Given Parameters ---

# Generation
s_per_generator_MVA = 100.0
num_generators = 3
power_factor = 0.9

# Loads
# There are 3 base loads of 50/3 MW each.
p_base_load_each_MW = 50.0 / 3.0
num_base_loads = 3
# The additional load is the active power one, 100 MW.
p_additional_load_MW = 100.0

# --- Calculations ---

# 1. Calculate Total Generated Active Power (P_gen)
# Total apparent power S = Number of generators * MVA per generator
s_gen_total_MVA = s_per_generator_MVA * num_generators
# Total active power P = Total S * Power Factor
p_gen_total_MW = s_gen_total_MVA * power_factor

# 2. Calculate Total Active Power Load (P_load)
# Total base load = Number of base loads * Power per base load
p_base_load_total_MW = p_base_load_each_MW * num_base_loads
# Total load = Total base load + Additional load
p_load_total_MW = p_base_load_total_MW + p_additional_load_MW

# 3. Calculate Total Resistive Power Loss (P_loss)
# Power Loss = Power Generated - Power Load
p_loss_MW = p_gen_total_MW - p_load_total_MW

# --- Output the results ---

print("Step 1: Calculate Total Generated Active Power (P_gen)")
print(f"Total Apparent Power (S_gen) = {s_per_generator_MVA:.1f} MVA/gen * {num_generators} gen = {s_gen_total_MVA:.1f} MVA")
print(f"P_gen = S_gen * Power Factor")
print(f"P_gen = {s_gen_total_MVA:.1f} MVA * {power_factor} = {p_gen_total_MW:.1f} MW")
print("-" * 40)

print("Step 2: Calculate Total Active Power Load (P_load)")
# Note: 50/3 MW * 3 loads = 50 MW exactly.
print(f"Total Base Load = (50/3 MW) * {num_base_loads} loads = {p_base_load_total_MW:.1f} MW")
print(f"Additional Load = {p_additional_load_MW:.1f} MW")
print(f"P_load = {p_base_load_total_MW:.1f} MW + {p_additional_load_MW:.1f} MW = {p_load_total_MW:.1f} MW")
print("-" * 40)

print("Step 3: Calculate Total Resistive Power Loss (P_loss)")
print("P_loss = P_gen - P_load")
print(f"P_loss = {p_gen_total_MW:.1f} MW - {p_load_total_MW:.1f} MW = {p_loss_MW:.1f} MW")
print("-" * 40)

print(f"The total power losses in the system are {p_loss_MW:.3f} MW.")