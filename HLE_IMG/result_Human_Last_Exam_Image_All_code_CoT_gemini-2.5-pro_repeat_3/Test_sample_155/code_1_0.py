import math

# Step 1: Define all the given values from the problem description and diagram.

# Load real power values (in MW)
p_load_1 = 1.89  # MW
p_load_2 = 1.75  # MW
p_load_3 = 1.7   # MW
p_load_4 = 1.9   # MW
p_load_5 = 2.4   # MW

# Generation apparent power values (in MVA)
s_hydro_gen = 2.0      # MVA
s_mini_hydro_gen = 1.2 # MVA
s_pv = 3.0           # MVA

# Power factors
pf_sources = 0.9  # for Hydro Gen and Mini Hydro Gen
pf_pv = 0.92      # for Photovoltaic system

# Losses
harmonic_loss = 0.05  # 5%
pv_transmission_loss = 0.04 # 4%

# Step 2: Calculate the total real power demand from all loads.

# Apply 5% harmonic loss to Load 2 and Load 4
p_load_2_eff = p_load_2 * (1 + harmonic_loss)
p_load_4_eff = p_load_4 * (1 + harmonic_loss)

# Sum of all loads
total_load_power = p_load_1 + p_load_2_eff + p_load_3 + p_load_4_eff + p_load_5

# Step 3: Calculate the total real power supplied by local generators.

# Real power from Hydro Generator (P = S * PF)
p_hydro_gen = s_hydro_gen * pf_sources

# Real power from Mini Hydro Generator (P = S * PF)
p_mini_hydro_gen = s_mini_hydro_gen * pf_sources

# Real power from Photovoltaic system, accounting for transmission loss
p_pv_generated = s_pv * pf_pv
p_pv_delivered = p_pv_generated * (1 - pv_transmission_loss)

# Sum of all generation
total_generation_power = p_hydro_gen + p_mini_hydro_gen + p_pv_delivered

# Step 4: Calculate the net real power demand on Bus 11.
net_real_power_demand = total_load_power - total_generation_power

# Step 5: Print the results and the final equation.

print("Calculation of Total Net Real Power Demand on Bus 11\n")

print("1. Calculating Total Load Power:")
print(f"  - Load 1: {p_load_1:.2f} MW")
print(f"  - Load 2 (with 5% harmonic loss): {p_load_2:.2f} MW * 1.05 = {p_load_2_eff:.4f} MW")
print(f"  - Load 3: {p_load_3:.2f} MW")
print(f"  - Load 4 (with 5% harmonic loss): {p_load_4:.2f} MW * 1.05 = {p_load_4_eff:.4f} MW")
print(f"  - Load 5: {p_load_5:.2f} MW")
print(f"  - Total Load Power = {p_load_1:.2f} + {p_load_2_eff:.4f} + {p_load_3:.2f} + {p_load_4_eff:.4f} + {p_load_5:.2f} = {total_load_power:.4f} MW\n")


print("2. Calculating Total Generation Power:")
print(f"  - Hydro Gen: {s_hydro_gen:.1f} MVA * {pf_sources} = {p_hydro_gen:.2f} MW")
print(f"  - Mini Hydro Gen: {s_mini_hydro_gen:.1f} MVA * {pf_sources} = {p_mini_hydro_gen:.2f} MW")
print(f"  - Photovoltaic (net): ({s_pv:.1f} MVA * {pf_pv}) * (1 - {pv_transmission_loss}) = {p_pv_delivered:.4f} MW")
print(f"  - Total Generation Power = {p_hydro_gen:.2f} + {p_mini_hydro_gen:.2f} + {p_pv_delivered:.4f} = {total_generation_power:.4f} MW\n")

print("3. Calculating Net Real Power Demand on Bus 11:")
print("   Net Demand = Total Load Power - Total Generation Power")
print(f"   Net Demand = ({p_load_1:.2f} + {p_load_2_eff:.4f} + {p_load_3:.1f} + {p_load_4_eff:.4f} + {p_load_5:.1f}) - ({p_hydro_gen:.1f} + {p_mini_hydro_gen:.2f} + {p_pv_delivered:.4f})")
print(f"   Net Demand = {total_load_power:.4f} MW - {total_generation_power:.4f} MW")
print(f"   Net Demand = {net_real_power_demand:.4f} MW")

# Final answer in the required format
#<<<4.2929>>>