import math

# --- Given Data ---

# Loads (Real Power in MW)
p_load1 = 1.89
p_load2 = 1.75
p_load3 = 1.7
p_load4 = 1.9
p_load5 = 2.4

# Generators (Apparent Power in MVA)
s_hydro_gen = 2.0
s_mini_hydro_gen = 1.2
s_pv_gen = 3.0

# Power Factors and Losses
pf_sources = 0.9  # For Hydro and Mini Hydro
pf_pv = 0.92
harmonic_loss_pct = 5.0
pv_transmission_loss_pct = 4.0

# --- Calculations ---

# 1. Calculate total load power demand
harmonic_loss_multiplier = 1 + harmonic_loss_pct / 100
p_load2_actual = p_load2 * harmonic_loss_multiplier
p_load4_actual = p_load4 * harmonic_loss_multiplier

total_load = p_load1 + p_load2_actual + p_load3 + p_load4_actual + p_load5

# 2. Calculate total generation power supplied
p_hydro_gen = s_hydro_gen * pf_sources
p_mini_hydro_gen = s_mini_hydro_gen * pf_sources

p_pv_generated = s_pv_gen * pf_pv
pv_transmission_loss_multiplier = 1 - pv_transmission_loss_pct / 100
p_pv_delivered = p_pv_generated * pv_transmission_loss_multiplier

total_generation = p_hydro_gen + p_mini_hydro_gen + p_pv_delivered

# 3. Calculate net real power demand on Bus 11
net_demand = total_load - total_generation

# --- Output the results and the equation ---
print("Step 1: Calculate Total Real Power Load")
print(f"Total Load = P_L1 + (P_L2 * (1 + {harmonic_loss_pct/100})) + P_L3 + (P_L4 * (1 + {harmonic_loss_pct/100})) + P_L5")
print(f"Total Load = {p_load1} MW + ({p_load2} MW * {harmonic_loss_multiplier}) + {p_load3} MW + ({p_load4} MW * {harmonic_loss_multiplier}) + {p_load5} MW")
print(f"Total Load = {p_load1} + {p_load2_actual:.4f} + {p_load3} + {p_load4_actual:.4f} + {p_load5} = {total_load:.4f} MW\n")

print("Step 2: Calculate Total Real Power Generation")
print(f"Total Generation = (S_hydro * PF_hydro) + (S_mini_hydro * PF_mini_hydro) + (S_pv * PF_pv * (1 - loss_pv))")
print(f"Total Generation = ({s_hydro_gen} MVA * {pf_sources}) + ({s_mini_hydro_gen} MVA * {pf_sources}) + ({s_pv_gen} MVA * {pf_pv} * {pv_transmission_loss_multiplier})")
print(f"Total Generation = {p_hydro_gen:.4f} MW + {p_mini_hydro_gen:.4f} MW + {p_pv_delivered:.4f} MW = {total_generation:.4f} MW\n")

print("Step 3: Calculate Net Real Power Demand on Bus 11")
print(f"Net Demand = Total Load - Total Generation")
print(f"Net Demand = {total_load:.4f} MW - {total_generation:.4f} MW")
print(f"Net Demand = {net_demand:.4f} MW")
