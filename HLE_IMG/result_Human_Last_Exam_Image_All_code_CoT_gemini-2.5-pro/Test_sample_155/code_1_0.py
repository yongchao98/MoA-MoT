import math

# Step 1: Define all the given values from the problem description and diagram.

# Load values (in MW)
load1_p = 1.89
load2_p_base = 1.75
load3_p = 1.7
load4_p_base = 1.9
load5_p = 2.4

# Generation values
hydro_gen_s = 2.0  # MVA
mini_hydro_gen_s = 1.2  # MVA
pv_s = 3.0  # MVA

# Power factors
general_source_pf = 0.9
pv_pf = 0.92

# Losses
harmonic_loss_pct = 0.05  # 5%
pv_transmission_loss_pct = 0.04 # 4%

# Step 2: Calculate the total real power load, accounting for losses.

# Calculate effective power for loads with harmonic loss
load2_p_eff = load2_p_base * (1 + harmonic_loss_pct)
load4_p_eff = load4_p_base * (1 + harmonic_loss_pct)

# Sum all loads to get total load power
total_load_p = load1_p + load2_p_eff + load3_p + load4_p_eff + load5_p

# Step 3: Calculate the total real power generation, accounting for losses.

# Calculate real power for each generator
hydro_gen_p = hydro_gen_s * general_source_pf
mini_hydro_gen_p = mini_hydro_gen_s * general_source_pf

# Calculate PV power generated and then the net power delivered after transmission loss
pv_p_generated = pv_s * pv_pf
pv_p_net = pv_p_generated * (1 - pv_transmission_loss_pct)

# Sum all generation to get total generation power
total_generation_p = hydro_gen_p + mini_hydro_gen_p + pv_p_net

# Step 4: Calculate the net real power demand on Bus 11.
net_demand_p = total_load_p - total_generation_p

# Step 5: Print the detailed calculation and the final answer.
print("Calculation of Total Net Real Power Demand on Bus 11")
print("------------------------------------------------------")

print("\n1. Total Real Power Load Calculation:")
print(f"  Load 1 = {load1_p:.2f} MW")
print(f"  Load 2 (effective) = {load2_p_base:.2f} MW * (1 + {harmonic_loss_pct}) = {load2_p_eff:.4f} MW")
print(f"  Load 3 = {load3_p:.2f} MW")
print(f"  Load 4 (effective) = {load4_p_base:.2f} MW * (1 + {harmonic_loss_pct}) = {load4_p_eff:.4f} MW")
print(f"  Load 5 = {load5_p:.2f} MW")
print(f"  Total Load = {load1_p:.2f} + {load2_p_eff:.4f} + {load3_p:.2f} + {load4_p_eff:.4f} + {load5_p:.2f} = {total_load_p:.4f} MW")

print("\n2. Total Real Power Generation Calculation:")
print(f"  Hydro Gen. = {hydro_gen_s:.1f} MVA * {general_source_pf} = {hydro_gen_p:.2f} MW")
print(f"  Mini Hydro Gen. = {mini_hydro_gen_s:.1f} MVA * {general_source_pf} = {mini_hydro_gen_p:.2f} MW")
print(f"  Photovoltaic (net) = ({pv_s:.1f} MVA * {pv_pf}) * (1 - {pv_transmission_loss_pct}) = {pv_p_net:.4f} MW")
print(f"  Total Generation = {hydro_gen_p:.2f} + {mini_hydro_gen_p:.2f} + {pv_p_net:.4f} = {total_generation_p:.4f} MW")

print("\n3. Net Real Power Demand Calculation:")
print("  Net Demand = Total Load - Total Generation")
print(f"  Net Demand = {total_load_p:.4f} MW - {total_generation_p:.4f} MW = {net_demand_p:.4f} MW")

print("\nFinal Equation:")
print(f"Total Net Real Power Demand = ({load1_p} + {load2_p_eff} + {load3_p} + {load4_p_eff} + {load5_p}) - ({hydro_gen_p} + {mini_hydro_gen_p} + {pv_p_net}) = {net_demand_p:.4f} MW")
