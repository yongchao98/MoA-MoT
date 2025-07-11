import math

# Step 1: Define initial values for loads and generators from the diagram and problem description.
p_load1 = 1.89  # MW
p_load2 = 1.75  # MW
p_load3 = 1.7   # MW
p_load4 = 1.9   # MW
p_load5 = 2.4   # MW

s_hydro_gen = 2.0      # MVA
s_mini_hydro_gen = 1.2 # MVA
s_pv = 3.0             # MVA

# Step 2: Define operational parameters.
pf_sources = 0.9  # Power factor for Hydro and Mini Hydro generators
pf_pv = 0.92      # Power factor for the Photovoltaic system
harmonic_loss_pct = 0.05  # 5% harmonic loss for Load 2 and Load 4
pv_transmission_loss_pct = 0.04  # 4% transmission loss for the PV system

# Step 3: Calculate the total real power demand from all loads.
# For loads with harmonic losses, the power drawn from the bus is P_load / (1 - loss).
p_load2_adj = p_load2 / (1 - harmonic_loss_pct)
p_load4_adj = p_load4 / (1 - harmonic_loss_pct)

total_load_power = p_load1 + p_load2_adj + p_load3 + p_load4_adj + p_load5

print("--- Calculating Total Load Power ---")
print(f"Total Load = P_Load1 + P_Load2_adj + P_Load3 + P_Load4_adj + P_Load5")
print(f"Total Load = {p_load1:.2f} MW + ({p_load2:.2f} / (1-{harmonic_loss_pct})) MW + {p_load3:.1f} MW + ({p_load4:.1f} / (1-{harmonic_loss_pct})) MW + {p_load5:.1f} MW")
print(f"Total Load = {p_load1:.2f} + {p_load2_adj:.4f} + {p_load3:.1f} + {p_load4_adj:.4f} + {p_load5:.1f} = {total_load_power:.4f} MW\n")


# Step 4: Calculate the total real power supplied by all generators to Bus 11.
# Real Power (P) = Apparent Power (S) * Power Factor (PF)
p_hydro_gen = s_hydro_gen * pf_sources
p_mini_hydro_gen = s_mini_hydro_gen * pf_sources

# For the PV system, calculate generated power first, then apply the transmission loss.
p_pv_generated = s_pv * pf_pv
p_pv_delivered = p_pv_generated * (1 - pv_transmission_loss_pct)

total_generation_power = p_hydro_gen + p_mini_hydro_gen + p_pv_delivered

print("--- Calculating Total Generation Power ---")
print(f"Total Generation = P_Hydro + P_MiniHydro + P_PV_delivered")
print(f"Total Generation = ({s_hydro_gen:.1f} MVA * {pf_sources}) + ({s_mini_hydro_gen:.1f} MVA * {pf_sources}) + (({s_pv:.1f} MVA * {pf_pv}) * (1-{pv_transmission_loss_pct}))")
print(f"Total Generation = {p_hydro_gen:.2f} + {p_mini_hydro_gen:.2f} + {p_pv_delivered:.4f} = {total_generation_power:.4f} MW\n")

# Step 5: Calculate the net real power demand on Bus 11.
net_demand = total_load_power - total_generation_power

print("--- Calculating Net Real Power Demand on Bus 11 ---")
print(f"Net Demand = Total Load Power - Total Generation Power")
print(f"Net Demand = {total_load_power:.4f} MW - {total_generation_power:.4f} MW")
print(f"Net Demand = {net_demand:.4f} MW")

# The final answer in the required format
# print(f"\n<<<{net_demand:.4f}>>>")