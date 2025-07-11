import math

# Define given values for loads
p_load1 = 1.89  # MW
p_load2_base = 1.75  # MW
p_load3 = 1.7   # MW
p_load4_base = 1.9   # MW
p_load5 = 2.4   # MW

# Define given values for generators
s_hydro = 2.0      # MVA, Hydro Gen.
s_mini_hydro = 1.2 # MVA, Mini Hydro Gen.
s_pv = 3.0         # MVA, Photovoltaic

# Define operational parameters
pf_sources = 0.9      # Power factor for most sources
pf_pv = 0.92          # Power factor for PV system
loss_harmonic_rate = 0.05  # 5%
loss_pv_transmission_rate = 0.04 # 4%

# --- Step 1: Calculate total load power ---
# Add harmonic loss for Load 2 and Load 4
p_load2_total = p_load2_base * (1 + loss_harmonic_rate)
p_load4_total = p_load4_base * (1 + loss_harmonic_rate)

# Sum all loads
total_load_power = p_load1 + p_load2_total + p_load3 + p_load4_total + p_load5

print("Step 1: Calculate Total Real Power Load")
print("---------------------------------------")
print(f"Load 2 with 5% harmonic loss: {p_load2_base:.2f} MW * (1 + {loss_harmonic_rate}) = {p_load2_total:.4f} MW")
print(f"Load 4 with 5% harmonic loss: {p_load4_base:.2f} MW * (1 + {loss_harmonic_rate}) = {p_load4_total:.4f} MW")
print("Total Load Power = Load 1 + Load 2 (adj) + Load 3 + Load 4 (adj) + Load 5")
print(f"Total Load Power = {p_load1:.2f} MW + {p_load2_total:.4f} MW + {p_load3:.1f} MW + {p_load4_total:.4f} MW + {p_load5:.1f} MW = {total_load_power:.4f} MW")
print("\n")


# --- Step 2: Calculate total generation power ---
# Calculate real power for each generator
p_hydro = s_hydro * pf_sources
p_mini_hydro = s_mini_hydro * pf_sources
p_pv_net = (s_pv * pf_pv) * (1 - loss_pv_transmission_rate)

# Sum all generation
total_generation_power = p_hydro + p_mini_hydro + p_pv_net

print("Step 2: Calculate Total Real Power Generation")
print("-------------------------------------------")
print(f"Hydro Gen. Power: {s_hydro:.1f} MVA * {pf_sources} = {p_hydro:.2f} MW")
print(f"Mini Hydro Gen. Power: {s_mini_hydro:.1f} MVA * {pf_sources} = {p_mini_hydro:.2f} MW")
print(f"PV System Net Power (after 4% loss): ({s_pv:.1f} MVA * {pf_pv}) * (1 - {loss_pv_transmission_rate}) = {p_pv_net:.4f} MW")
print("Total Generation Power = Hydro Gen + Mini Hydro Gen + PV (net)")
print(f"Total Generation Power = {p_hydro:.2f} MW + {p_mini_hydro:.2f} MW + {p_pv_net:.4f} MW = {total_generation_power:.4f} MW")
print("\n")


# --- Step 3: Calculate net real power demand ---
net_demand = total_load_power - total_generation_power

print("Step 3: Calculate Net Real Power Demand on Bus 11")
print("-------------------------------------------------")
print("Net Demand = Total Load Power - Total Generation Power")
print(f"Net Demand = {total_load_power:.4f} MW - {total_generation_power:.4f} MW = {net_demand:.4f} MW")

# Final answer in the specified format
# print(f"<<<{net_demand:.4f}>>>")