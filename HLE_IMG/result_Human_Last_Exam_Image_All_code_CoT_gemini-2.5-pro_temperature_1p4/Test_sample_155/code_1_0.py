import math

# --- Given Data ---

# Loads (Real Power in MW)
p_load1 = 1.89
p_load2_base = 1.75
p_load3 = 1.7
p_load4_base = 1.9
p_load5 = 2.4

# Generators (Apparent Power in MVA)
s_hydro_gen = 2.0
s_mini_hydro_gen = 1.2
s_pv_gen = 3.0

# System Parameters
pf_sources = 0.9  # For Hydro and Mini Hydro
pf_pv = 0.92
harmonic_loss_pct = 5.0
transmission_loss_pv_pct = 4.0

# --- Step-by-step Calculation ---

print("Step 1: Calculate the effective real power demand from each load.")
# Calculate power demand for loads with harmonic losses
p_load2_actual = p_load2_base * (1 + harmonic_loss_pct / 100)
p_load4_actual = p_load4_base * (1 + harmonic_loss_pct / 100)
print(f"Load 1 Demand: {p_load1:.4f} MW")
print(f"Load 2 Demand (with {harmonic_loss_pct}% harmonic loss): {p_load2_base:.2f} MW * 1.05 = {p_load2_actual:.4f} MW")
print(f"Load 3 Demand: {p_load3:.4f} MW")
print(f"Load 4 Demand (with {harmonic_loss_pct}% harmonic loss): {p_load4_base:.2f} MW * 1.05 = {p_load4_actual:.4f} MW")
print(f"Load 5 Demand: {p_load5:.4f} MW")

# Calculate total load
total_load = p_load1 + p_load2_actual + p_load3 + p_load4_actual + p_load5
print(f"Total Real Power Load = {p_load1:.2f} + {p_load2_actual:.4f} + {p_load3:.2f} + {p_load4_actual:.4f} + {p_load5:.2f} = {total_load:.4f} MW\n")


print("Step 2: Calculate the real power supplied by each generator to Bus 11.")
# Calculate real power from generators
p_hydro_gen = s_hydro_gen * pf_sources
p_mini_hydro_gen = s_mini_hydro_gen * pf_sources

# Calculate PV power delivered after losses
p_pv_generated = s_pv_gen * pf_pv
p_pv_delivered = p_pv_generated * (1 - transmission_loss_pv_pct / 100)

print(f"Hydro Gen. Supply: {s_hydro_gen:.1f} MVA * {pf_sources} PF = {p_hydro_gen:.4f} MW")
print(f"Mini Hydro Gen. Supply: {s_mini_hydro_gen:.1f} MVA * {pf_sources} PF = {p_mini_hydro_gen:.4f} MW")
print(f"Photovoltaic Supply (after {transmission_loss_pv_pct}% loss): ({s_pv_gen:.1f} MVA * {pf_pv} PF) * (1 - {transmission_loss_pv_pct/100}) = {p_pv_delivered:.4f} MW")

# Calculate total generation
total_generation = p_hydro_gen + p_mini_hydro_gen + p_pv_delivered
print(f"Total Real Power Generation = {p_hydro_gen:.4f} + {p_mini_hydro_gen:.4f} + {p_pv_delivered:.4f} = {total_generation:.4f} MW\n")


print("Step 3: Calculate the net real power demand on Bus 11.")
# Calculate net demand
net_demand = total_load - total_generation

# Print the final equation with all numbers
print("Final Equation:")
print(f"Net Demand = (P_L1 + P_L2_actual + P_L3 + P_L4_actual + P_L5) - (P_Hydro + P_MiniHydro + P_PV_delivered)")
print(f"Net Demand = ({p_load1:.2f} + {p_load2_actual:.4f} + {p_load3:.2f} + {p_load4_actual:.4f} + {p_load5:.2f}) - ({p_hydro_gen:.4f} + {p_mini_hydro_gen:.4f} + {p_pv_delivered:.4f})")
print(f"Net Demand = {total_load:.4f} MW - {total_generation:.4f} MW = {net_demand:.4f} MW\n")

print("The total net real power demand on Bus 11 is:")
print(f"{net_demand:.4f} MW")
print("<<<4.2929>>>")