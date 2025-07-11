import math

# --- Step 1: Define initial values from the problem description ---

# Load values (in MW)
p_load1 = 1.89
p_load2_base = 1.75
p_load3 = 1.7
p_load4_base = 1.9
p_load5 = 2.4

# Generation values
s_hydro_gen = 2.0  # MVA
s_mini_hydro_gen = 1.2  # MVA
s_pv = 3.0  # MVA

# Power factors
pf_sources = 0.9  # For Hydro and Mini Hydro generators
pf_pv = 0.92      # For Photovoltaic system

# Losses
harmonic_loss_pct = 0.05  # 5% for Load 2 and Load 4
pv_transmission_loss_pct = 0.04  # 4% for PV system

# --- Step 2: Calculate total real power load ---

# Adjust loads for harmonic losses
p_load2_actual = p_load2_base * (1 + harmonic_loss_pct)
p_load4_actual = p_load4_base * (1 + harmonic_loss_pct)

# Sum all loads
total_load_p = p_load1 + p_load2_actual + p_load3 + p_load4_actual + p_load5

print("--- Calculating Total Real Power Load ---")
print(f"Load 1 Power: {p_load1:.2f} MW")
print(f"Load 2 Power (including 5% harmonic loss): {p_load2_base:.2f} * {1 + harmonic_loss_pct} = {p_load2_actual:.4f} MW")
print(f"Load 3 Power: {p_load3:.2f} MW")
print(f"Load 4 Power (including 5% harmonic loss): {p_load4_base:.2f} * {1 + harmonic_loss_pct} = {p_load4_actual:.4f} MW")
print(f"Load 5 Power: {p_load5:.2f} MW")
print(f"Total Real Power Load = {p_load1:.2f} + {p_load2_actual:.4f} + {p_load3:.2f} + {p_load4_actual:.4f} + {p_load5:.2f} = {total_load_p:.4f} MW\n")


# --- Step 3: Calculate total real power generation ---

# Calculate real power for each generator
p_hydro_gen = s_hydro_gen * pf_sources
p_mini_hydro_gen = s_mini_hydro_gen * pf_sources
p_pv_generated = s_pv * pf_pv

# Adjust PV generation for transmission loss
p_pv_net = p_pv_generated * (1 - pv_transmission_loss_pct)

# Sum all generations
total_generation_p = p_hydro_gen + p_mini_hydro_gen + p_pv_net

print("--- Calculating Total Real Power Generation ---")
print(f"Hydro Gen. Power: {s_hydro_gen:.1f} MVA * {pf_sources} PF = {p_hydro_gen:.2f} MW")
print(f"Mini Hydro Gen. Power: {s_mini_hydro_gen:.1f} MVA * {pf_sources} PF = {p_mini_hydro_gen:.2f} MW")
print(f"Net Photovoltaic Power (after 4% loss): ({s_pv:.1f} MVA * {pf_pv} PF) * (1 - {pv_transmission_loss_pct}) = {p_pv_net:.4f} MW")
print(f"Total Real Power Generation = {p_hydro_gen:.2f} + {p_mini_hydro_gen:.2f} + {p_pv_net:.4f} = {total_generation_p:.4f} MW\n")

# --- Step 4: Calculate net real power demand on Bus 11 ---

net_demand_p = total_load_p - total_generation_p

print("--- Calculating Net Real Power Demand on Bus 11 ---")
print("Net Demand = Total Load - Total Generation")
print(f"Net Demand = ({p_load1:.2f} + {p_load2_actual:.4f} + {p_load3:.2f} + {p_load4_actual:.4f} + {p_load5:.2f}) - ({p_hydro_gen:.2f} + {p_mini_hydro_gen:.2f} + {p_pv_net:.4f})")
print(f"Net Demand = {total_load_p:.4f} MW - {total_generation_p:.4f} MW")
print(f"Final Net Real Power Demand on Bus 11 = {net_demand_p:.4f} MW")

# Final answer in the required format
final_answer = net_demand_p
print(f"\n<<<{final_answer:.4f}>>>")