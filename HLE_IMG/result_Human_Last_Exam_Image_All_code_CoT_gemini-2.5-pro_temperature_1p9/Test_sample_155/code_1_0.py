import math

# --- 1. Define given values from the microgrid diagram ---

# Load real power values (in MW)
load1_p_rated = 1.89
load2_p_rated = 1.75
load3_p_rated = 1.7
load4_p_rated = 1.9
load5_p_rated = 2.4

# Generator apparent power values (in MVA)
hydro_gen_s = 2.0
mini_hydro_gen_s = 1.2
pv_s = 3.0

# Power factors (PF)
# For "Load 3 and all sources"
general_pf = 0.9
# Specific PF for the PV system
pv_pf = 0.92

# Losses
harmonic_loss_percent = 5.0
pv_transmission_loss_percent = 4.0

# --- 2. Calculate the total load demand ---

print("### Calculating Total Real Power Load ###\n")

# Calculate actual power for loads with harmonic losses
harmonic_multiplier = 1 + (harmonic_loss_percent / 100)
load2_p_actual = load2_p_rated * harmonic_multiplier
load4_p_actual = load4_p_rated * harmonic_multiplier

print(f"Load 2 Power with 5% harmonic loss: {load2_p_rated} MW * {harmonic_multiplier} = {load2_p_actual:.4f} MW")
print(f"Load 4 Power with 5% harmonic loss: {load4_p_rated} MW * {harmonic_multiplier} = {load4_p_actual:.4f} MW\n")

# Sum all loads
total_load_p = load1_p_rated + load2_p_actual + load3_p_rated + load4_p_actual + load5_p_rated
print("Total Load Demand Calculation:")
print(f"Total Load = Load 1 + Load 2 (actual) + Load 3 + Load 4 (actual) + Load 5")
print(f"Total Load = {load1_p_rated} MW + {load2_p_actual:.4f} MW + {load3_p_rated} MW + {load4_p_actual:.4f} MW + {load5_p_rated} MW")
print(f"Total Real Power Load = {total_load_p:.4f} MW\n")

# --- 3. Calculate the total power generation ---

print("### Calculating Total Real Power Generation ###\n")

# Calculate real power for each generator
hydro_gen_p = hydro_gen_s * general_pf
mini_hydro_gen_p = mini_hydro_gen_s * general_pf
pv_p_initial = pv_s * pv_pf

# Apply transmission loss to PV generation
pv_loss_multiplier = 1 - (pv_transmission_loss_percent / 100)
pv_p_actual = pv_p_initial * pv_loss_multiplier

print(f"Hydro Gen. Power: {hydro_gen_s} MVA * {general_pf} PF = {hydro_gen_p:.4f} MW")
print(f"Mini Hydro Gen. Power: {mini_hydro_gen_s} MVA * {general_pf} PF = {mini_hydro_gen_p:.4f} MW")
print(f"Photovoltaic Power (before loss): {pv_s} MVA * {pv_pf} PF = {pv_p_initial:.4f} MW")
print(f"Photovoltaic Power (after 4% loss): {pv_p_initial:.4f} MW * {pv_loss_multiplier} = {pv_p_actual:.4f} MW\n")

# Sum all generated power
total_generation_p = hydro_gen_p + mini_hydro_gen_p + pv_p_actual
print("Total Generation Calculation:")
print(f"Total Generation = Hydro Gen. + Mini Hydro Gen. + Photovoltaic (actual)")
print(f"Total Generation = {hydro_gen_p:.4f} MW + {mini_hydro_gen_p:.4f} MW + {pv_p_actual:.4f} MW")
print(f"Total Real Power Generation = {total_generation_p:.4f} MW\n")


# --- 4. Calculate the net demand on Bus 11 ---

print("### Calculating Net Real Power Demand on Bus 11 ###\n")

net_demand_p = total_load_p - total_generation_p

print("Final Equation:")
print(f"Net Demand = Total Load - Total Generation")
print(f"Net Demand = {total_load_p:.4f} MW - {total_generation_p:.4f} MW")
print(f"Total Net Real Power Demand on Bus 11 = {net_demand_p:.4f} MW")