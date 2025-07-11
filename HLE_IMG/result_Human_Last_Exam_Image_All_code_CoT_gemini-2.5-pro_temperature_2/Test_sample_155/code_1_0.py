import math

# Step 1: Define all the given values from the problem description
# --- Loads (in MW) ---
p_load1 = 1.89  # MW
p_load2_base = 1.75  # MW
p_load3 = 1.7  # MW
p_load4_base = 1.9  # MW
p_load5 = 2.4  # MW
harmonic_loss_pct = 5.0 # %

# --- Generators (Apparent Power in MVA, Power Factor) ---
s_hydro_gen = 2.0  # MVA
pf_sources = 0.9  # Lagging for Hydro and Mini Hydro

s_mini_hydro_gen = 1.2 # MVA

s_pv_gen = 3.0  # MVA
pf_pv = 0.92  # Lagging
pv_transmission_loss_pct = 4.0 # %

# Step 2: Calculate the total real power load
# Account for harmonic losses in Load 2 and Load 4
p_load2_actual = p_load2_base * (1 + harmonic_loss_pct / 100.0)
p_load4_actual = p_load4_base * (1 + harmonic_loss_pct / 100.0)

# Sum all loads
total_load_power = p_load1 + p_load2_actual + p_load3 + p_load4_actual + p_load5

# Step 3: Calculate the total real power generation delivered to Bus 11
# Calculate real power for each generator
p_hydro_gen = s_hydro_gen * pf_sources
p_mini_hydro_gen = s_mini_hydro_gen * pf_sources

# Calculate PV power generated and then delivered power after losses
p_pv_generated = s_pv_gen * pf_pv
p_pv_delivered = p_pv_generated * (1 - pv_transmission_loss_pct / 100.0)

# Sum all generated powers
total_generation_power = p_hydro_gen + p_mini_hydro_gen + p_pv_delivered

# Step 4: Calculate the net real power demand on Bus 11
net_power_demand = total_load_power - total_generation_power

# --- Output the results step-by-step ---
print("--- Calculation of Net Real Power Demand on Bus 11 ---")
print("\n1. Calculating Total Load Power:")
print(f"  - Load 1: {p_load1:.2f} MW")
print(f"  - Load 2: {p_load2_base:.2f} MW + {harmonic_loss_pct}% harmonic loss = {p_load2_actual:.4f} MW")
print(f"  - Load 3: {p_load3:.2f} MW")
print(f"  - Load 4: {p_load4_base:.2f} MW + {harmonic_loss_pct}% harmonic loss = {p_load4_actual:.4f} MW")
print(f"  - Load 5: {p_load5:.2f} MW")
print(f"  - Total Load Power = {p_load1:.2f} + {p_load2_actual:.4f} + {p_load3:.2f} + {p_load4_actual:.4f} + {p_load5:.2f} = {total_load_power:.4f} MW")

print("\n2. Calculating Total Generation Power Delivered to Bus 11:")
print(f"  - Hydro Gen: {s_hydro_gen:.1f} MVA * {pf_sources:.2f} PF = {p_hydro_gen:.2f} MW")
print(f"  - Mini Hydro Gen: {s_mini_hydro_gen:.1f} MVA * {pf_sources:.2f} PF = {p_mini_hydro_gen:.2f} MW")
print(f"  - Photovoltaic: ({s_pv_gen:.1f} MVA * {pf_pv:.2f} PF) * (1 - {pv_transmission_loss_pct/100:.2f} loss) = {p_pv_delivered:.4f} MW")
print(f"  - Total Generation Power = {p_hydro_gen:.2f} + {p_mini_hydro_gen:.2f} + {p_pv_delivered:.4f} = {total_generation_power:.4f} MW")

print("\n3. Final Calculation of Net Demand:")
print(f"Net Real Power Demand = Total Load - Total Generation")
# The final equation with each number printed out
print(f"Net Demand = ({p_load1:.2f} + {p_load2_actual:.4f} + {p_load3:.2f} + {p_load4_actual:.4f} + {p_load5:.2f}) - ({p_hydro_gen:.2f} + {p_mini_hydro_gen:.2f} + {p_pv_delivered:.4f})")
print(f"Net Demand = {total_load_power:.4f} MW - {total_generation_power:.4f} MW")
print(f"\nFinal Net Real Power Demand on Bus 11: {net_power_demand:.4f} MW")
