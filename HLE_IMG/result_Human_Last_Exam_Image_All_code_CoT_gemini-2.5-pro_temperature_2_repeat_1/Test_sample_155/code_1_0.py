import math

# Step 1: Define initial values from the problem description and diagram

# Load real power in MW
load1_p = 1.89
load2_p = 1.75
load3_p = 1.7
load4_p = 1.9
load5_p = 2.4

# Generator apparent power in MVA
hydro_gen_s = 2.0
mini_hydro_gen_s = 1.2
pv_s = 3.0
bess_s = 3.0 # Assumed from the transformer rating on its line

# Power factors
general_source_pf = 0.9
pv_pf = 0.92

# Losses in percentage
harmonic_loss_pct = 5.0
pv_transmission_loss_pct = 4.0

# Step 2: Calculate the real power drawn by each load considering losses

# For Load 2 and 4, there is a 5% harmonic power loss
p_load2_actual = load2_p * (1 + harmonic_loss_pct / 100)
p_load4_actual = load4_p * (1 + harmonic_loss_pct / 100)

# Total real power demand from all loads
total_p_demand = load1_p + p_load2_actual + load3_p + p_load4_actual + load5_p

print("Calculation of Total Load Demand (P_demand):")
print(f"P_demand = P_L1 + P_L2_with_loss + P_L3 + P_L4_with_loss + P_L5")
print(f"P_demand = {load1_p:.2f} MW + ({load2_p:.2f} MW * {1 + harmonic_loss_pct / 100:.2f}) + {load3_p:.2f} MW + ({load4_p:.2f} MW * {1 + harmonic_loss_pct / 100:.2f}) + {load5_p:.2f} MW")
print(f"P_demand = {load1_p:.2f} MW + {p_load2_actual:.4f} MW + {load3_p:.2f} MW + {p_load4_actual:.4f} MW + {load5_p:.2f} MW")
print(f"Total Load Demand = {total_p_demand:.4f} MW\n")


# Step 3: Calculate the real power supplied by each generator reaching Bus 11

# Real power from Hydro Generator
p_hydro = hydro_gen_s * general_source_pf

# Real power from Mini Hydro Generator
p_mini_hydro = mini_hydro_gen_s * general_source_pf

# Real power from BESS
p_bess = bess_s * general_source_pf

# Real power from Photovoltaic system, accounting for transmission loss
p_pv_generated = pv_s * pv_pf
p_pv_delivered = p_pv_generated * (1 - pv_transmission_loss_pct / 100)

# Total real power supply from all generators
total_p_supply = p_hydro + p_mini_hydro + p_bess + p_pv_delivered

print("Calculation of Total Generation Supply (P_supply):")
print(f"P_supply = P_hydro + P_mini_hydro + P_BESS + P_PV_delivered")
print(f"P_supply = ({hydro_gen_s:.1f} MVA * {general_source_pf:.1f}) + ({mini_hydro_gen_s:.1f} MVA * {general_source_pf:.1f}) + ({bess_s:.1f} MVA * {general_source_pf:.1f}) + (({pv_s:.1f} MVA * {pv_pf:.2f}) * (1 - {pv_transmission_loss_pct/100:.2f}))")
print(f"P_supply = {p_hydro:.2f} MW + {p_mini_hydro:.2f} MW + {p_bess:.2f} MW + {p_pv_delivered:.4f} MW")
print(f"Total Generation Supply = {total_p_supply:.4f} MW\n")


# Step 4: Calculate the total net real power demand on Bus 11
net_p_demand_bus11 = total_p_demand - total_p_supply

print("Calculation of Net Real Power Demand on Bus 11:")
print(f"Net Demand = Total Load Demand - Total Generation Supply")
print(f"Net Demand = {total_p_demand:.4f} MW - {total_p_supply:.4f} MW")
print(f"Net Real Power Demand on Bus 11 = {net_p_demand_bus11:.4f} MW")
print(f">>>{net_p_demand_bus11:.4f}<<<")