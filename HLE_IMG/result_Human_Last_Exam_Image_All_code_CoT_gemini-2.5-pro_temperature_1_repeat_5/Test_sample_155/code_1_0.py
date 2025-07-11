import math

# --- Given Parameters ---

# Loads (in MW)
p_load1 = 1.89
p_load2_rated = 1.75
p_load3 = 1.7
p_load4_rated = 1.9
p_load5 = 2.4

# Generators (Apparent Power in MVA)
s_hydro = 2.0
s_mini_hydro = 1.2
s_pv = 3.0

# Power Factors
pf_sources = 0.9  # For Hydro, Mini Hydro, and Load 3
pf_pv = 0.92

# Losses
harmonic_loss_percentage = 5.0
pv_transmission_loss_percentage = 4.0

# --- Calculation ---

# Step 1: Calculate Total Real Power Load
# Adjust for harmonic losses where applicable
harmonic_loss_factor = harmonic_loss_percentage / 100.0
p_load2_actual = p_load2_rated / (1 - harmonic_loss_factor)
p_load4_actual = p_load4_rated / (1 - harmonic_loss_factor)

total_load = p_load1 + p_load2_actual + p_load3 + p_load4_actual + p_load5

# Step 2: Calculate Total Real Power Generation
# Calculate real power for each source
p_hydro = s_hydro * pf_sources
p_mini_hydro = s_mini_hydro * pf_sources
p_pv_generated = s_pv * pf_pv

# Adjust PV generation for transmission loss
pv_transmission_loss_factor = pv_transmission_loss_percentage / 100.0
p_pv_net = p_pv_generated * (1 - pv_transmission_loss_factor)

total_generation = p_hydro + p_mini_hydro + p_pv_net

# Step 3: Calculate Net Real Power Demand
net_demand = total_load - total_generation

# --- Output the results ---

print("Calculation of the Net Real Power Demand on Bus 11:\n")

print("1. Total Load Demand Calculation:")
print(f"   - Load 1: {p_load1:.3f} MW")
print(f"   - Load 2 (adj. for {harmonic_loss_percentage}% loss): {p_load2_rated:.2f} / (1 - {harmonic_loss_factor}) = {p_load2_actual:.3f} MW")
print(f"   - Load 3: {p_load3:.3f} MW")
print(f"   - Load 4 (adj. for {harmonic_loss_percentage}% loss): {p_load4_rated:.2f} / (1 - {harmonic_loss_factor}) = {p_load4_actual:.3f} MW")
print(f"   - Load 5: {p_load5:.3f} MW")
print(f"   Total Load = {p_load1:.3f} + {p_load2_actual:.3f} + {p_load3:.3f} + {p_load4_actual:.3f} + {p_load5:.3f} = {total_load:.3f} MW\n")

print("2. Total Generation Supply Calculation:")
print(f"   - Hydro Gen: {s_hydro:.1f} MVA * {pf_sources} = {p_hydro:.3f} MW")
print(f"   - Mini Hydro Gen: {s_mini_hydro:.1f} MVA * {pf_sources} = {p_mini_hydro:.3f} MW")
print(f"   - PV System (net of {pv_transmission_loss_percentage}% loss): ({s_pv:.1f} MVA * {pf_pv}) * (1 - {pv_transmission_loss_factor}) = {p_pv_net:.3f} MW")
print(f"   Total Generation = {p_hydro:.3f} + {p_mini_hydro:.3f} + {p_pv_net:.3f} = {total_generation:.3f} MW\n")

print("3. Final Net Demand Calculation:")
print("   Net Demand = Total Load - Total Generation")
print(f"   Net Demand = {total_load:.3f} MW - {total_generation:.3f} MW = {net_demand:.3f} MW")
print("\nFinal Equation:")
print(f"({p_load1:.2f} + {p_load2_actual:.3f} + {p_load3:.2f} + {p_load4_actual:.3f} + {p_load5:.2f}) - ({p_hydro:.3f} + {p_mini_hydro:.3f} + {p_pv_net:.3f}) = {net_demand:.3f} MW")

# Final Answer
final_answer = net_demand
# print(f'<<<{final_answer:.4f}>>>')