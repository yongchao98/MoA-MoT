import math

# --- Step 1: Define initial values from the problem statement and diagram ---

# Loads (in MW)
P_L1 = 1.89
P_L2_nominal = 1.75
P_L3 = 1.7
P_L4_nominal = 1.9
P_L5 = 2.4

# Generators (Apparent Power in MVA, Power Factor)
S_PV = 3.0
PF_PV = 0.92
S_Hydro = 2.0
S_MiniHydro = 1.2
PF_sources = 0.9  # for Hydro and Mini Hydro

# Losses and other factors
harmonic_loss_factor = 0.05
pv_transmission_loss_factor = 0.04

# --- Step 2: Calculate Total Real Power Load ---

# Adjust for harmonic losses for Load 2 and Load 4
P_L2_actual = P_L2_nominal * (1 + harmonic_loss_factor)
P_L4_actual = P_L4_nominal * (1 + harmonic_loss_factor)

# Sum all loads to get the total load demand
P_load_total = P_L1 + P_L2_actual + P_L3 + P_L4_actual + P_L5

print("--- Calculating Total Real Power Load ---")
print(f"Load 1 Real Power: {P_L1} MW")
print(f"Load 2 Real Power (including 5% harmonic loss): {P_L2_nominal} * (1 + {harmonic_loss_factor}) = {P_L2_actual:.4f} MW")
print(f"Load 3 Real Power: {P_L3} MW")
print(f"Load 4 Real Power (including 5% harmonic loss): {P_L4_nominal} * (1 + {harmonic_loss_factor}) = {P_L4_actual:.4f} MW")
print(f"Load 5 Real Power: {P_L5} MW")
print(f"Total Real Power Load = {P_load_total:.4f} MW\n")


# --- Step 3: Calculate Total Real Power Generation ---

# Calculate generated power at the source for each generator
P_PV_generated = S_PV * PF_PV
P_Hydro_generated = S_Hydro * PF_sources
P_MiniHydro_generated = S_MiniHydro * PF_sources

# Account for transmission loss for the PV system
P_PV_delivered = P_PV_generated * (1 - pv_transmission_loss_factor)

# Sum all delivered generation to get the total generation
P_gen_total = P_PV_delivered + P_Hydro_generated + P_MiniHydro_generated

print("--- Calculating Total Real Power Generation ---")
print(f"PV Power Generated at Source: {S_PV} MVA * {PF_PV} pf = {P_PV_generated:.4f} MW")
print(f"PV Power Delivered to Bus 11 (after 4% loss): {P_PV_generated:.4f} * (1 - {pv_transmission_loss_factor}) = {P_PV_delivered:.4f} MW")
print(f"Hydro Gen Power: {S_Hydro} MVA * {PF_sources} pf = {P_Hydro_generated:.4f} MW")
print(f"Mini Hydro Gen Power: {S_MiniHydro} MVA * {PF_sources} pf = {P_MiniHydro_generated:.4f} MW")
print(f"Total Real Power Generation = {P_gen_total:.4f} MW\n")


# --- Step 4: Calculate Net Real Power Demand on Bus 11 ---

# Net Demand = Total Load - Total Generation
net_demand = P_load_total - P_gen_total

print("--- Calculating Net Real Power Demand on Bus 11 ---")
print("Net Demand = Total Load - Total Generation")
print(f"The final equation with all values is:")
print(f"Net Demand = ({P_L1} + {P_L2_actual:.4f} + {P_L3} + {P_L4_actual:.4f} + {P_L5}) - ({P_PV_delivered:.4f} + {P_Hydro_generated:.4f} + {P_MiniHydro_generated:.4f})")
print(f"Net Demand = ({P_load_total:.4f}) - ({P_gen_total:.4f}) = {net_demand:.4f} MW")
print(f"\nFinal Answer: The total net real power demand on Bus 11 is {net_demand:.4f} MW.")