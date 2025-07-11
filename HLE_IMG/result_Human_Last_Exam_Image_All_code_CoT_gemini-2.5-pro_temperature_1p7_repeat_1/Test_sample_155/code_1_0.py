import math

# --- 1. Define given values from the problem description and diagram ---

# Loads (Real Power in MW)
P_load1 = 1.89
P_load2_base = 1.75
P_load3 = 1.7
P_load4_base = 1.9
P_load5 = 2.4

# Generators (Apparent Power in MVA and Power Factor)
S_hydro_gen = 2.0
PF_hydro_gen = 0.9  # As per "all sources operate with a power factor of 0.9 lagging"

S_mini_hydro_gen = 1.2
PF_mini_hydro_gen = 0.9  # As per "all sources operate with a power factor of 0.9 lagging"

S_pv = 3.0
PF_pv = 0.92

# Losses
harmonic_loss_percent = 5.0
pv_transmission_loss_percent = 4.0

# --- 2. Calculate Total Load Demand with Losses ---

print("Step 1: Calculate the total effective real power demand from loads.")

# Apply 5% harmonic loss to Load 2 and Load 4
harmonic_loss_multiplier = 1 + harmonic_loss_percent / 100.0
P_load2_effective = P_load2_base * harmonic_loss_multiplier
P_load4_effective = P_load4_base * harmonic_loss_multiplier

# Sum all effective load powers
total_load_power = P_load1 + P_load2_effective + P_load3 + P_load4_effective + P_load5

print(f"Total Load Power = Load 1 + Load 2 (eff) + Load 3 + Load 4 (eff) + Load 5")
print(f"Total Load Power = {P_load1:.2f} MW + ({P_load2_base:.2f} * {harmonic_loss_multiplier:.2f}) MW + {P_load3:.2f} MW + ({P_load4_base:.2f} * {harmonic_loss_multiplier:.2f}) MW + {P_load5:.2f} MW")
print(f"Total Load Power = {P_load1:.2f} MW + {P_load2_effective:.4f} MW + {P_load3:.2f} MW + {P_load4_effective:.4f} MW + {P_load5:.2f} MW = {total_load_power:.4f} MW\n")


# --- 3. Calculate Total Power Generation with Losses ---

print("Step 2: Calculate the total net real power supplied by generators.")

# Calculate real power for each generator
P_hydro_gen = S_hydro_gen * PF_hydro_gen
P_mini_hydro_gen = S_mini_hydro_gen * PF_mini_hydro_gen
P_pv_generated = S_pv * PF_pv

# Apply 4% transmission loss to the PV system
pv_loss_multiplier = 1 - pv_transmission_loss_percent / 100.0
P_pv_net = P_pv_generated * pv_loss_multiplier

# Sum all net generated power
total_generation_power = P_hydro_gen + P_mini_hydro_gen + P_pv_net

print(f"Total Generation Power = Hydro Gen + Mini Hydro Gen + Photovoltaic (net)")
print(f"Total Generation Power = ({S_hydro_gen:.1f} MVA * {PF_hydro_gen:.1f}) + ({S_mini_hydro_gen:.1f} MVA * {PF_mini_hydro_gen:.1f}) + ({S_pv:.1f} MVA * {PF_pv:.2f} * {pv_loss_multiplier:.2f})")
print(f"Total Generation Power = {P_hydro_gen:.2f} MW + {P_mini_hydro_gen:.2f} MW + {P_pv_net:.4f} MW = {total_generation_power:.4f} MW\n")


# --- 4. Calculate Net Real Power Demand on Bus 11 ---
net_demand_on_bus11 = total_load_power - total_generation_power

print("Step 3: Calculate the net real power demand on Bus 11.")
print("Net Demand = Total Load Power - Total Generation Power")
print(f"Net Demand = {total_load_power:.4f} MW - {total_generation_power:.4f} MW")
print(f"Final Answer: The total net real power demand is {net_demand_on_bus11:.4f} MW.\n")

print("The complete equation is:")
print(f"({P_load1} + {P_load2_base}*1.05 + {P_load3} + {P_load4_base}*1.05 + {P_load5}) - ({S_hydro_gen}*{PF_hydro_gen} + {S_mini_hydro_gen}*{PF_mini_hydro_gen} + {S_pv}*{PF_pv}*(1 - {pv_transmission_loss_percent/100})) = {net_demand_on_bus11:.4f}")
