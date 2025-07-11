import math

# --- Given Data ---
# Loads (P in MW)
P_load1 = 1.89
P_load2 = 1.75
P_load3 = 1.7
P_load4 = 1.9
P_load5 = 2.4

# Generators (S in MVA)
S_hydro_gen = 2.0
S_mini_hydro_gen = 1.2
S_pv = 3.0

# Power Factors
pf_sources = 0.9  # For Hydro, Mini Hydro, and Load 3
pf_pv = 0.92

# Losses (%)
harmonic_loss_pct = 5.0
pv_transmission_loss_pct = 4.0

# --- Calculation ---

# Step 1: Calculate total real power load, accounting for losses
print("Step 1: Calculating Total Real Power Load\n")

# Load 1 has no special conditions
P_load1_actual = P_load1
print(f"Load 1 Real Power = {P_load1_actual:.4f} MW")

# Load 2 has a 5% harmonic loss
harmonic_loss_factor_2 = 1 + harmonic_loss_pct / 100
P_load2_actual = P_load2 * harmonic_loss_factor_2
print(f"Load 2 Real Power (with {harmonic_loss_pct}% harmonic loss) = {P_load2:.2f} MW * {harmonic_loss_factor_2:.2f} = {P_load2_actual:.4f} MW")

# Load 3 power is given in MW, PF is irrelevant for real power
P_load3_actual = P_load3
print(f"Load 3 Real Power = {P_load3_actual:.4f} MW")

# Load 4 has a 5% harmonic loss
harmonic_loss_factor_4 = 1 + harmonic_loss_pct / 100
P_load4_actual = P_load4 * harmonic_loss_factor_4
print(f"Load 4 Real Power (with {harmonic_loss_pct}% harmonic loss) = {P_load4:.2f} MW * {harmonic_loss_factor_4:.2f} = {P_load4_actual:.4f} MW")

# Load 5 reactive power is irrelevant for real power calculation
P_load5_actual = P_load5
print(f"Load 5 Real Power = {P_load5_actual:.4f} MW")

# Total Load Power
total_load_power = P_load1_actual + P_load2_actual + P_load3_actual + P_load4_actual + P_load5_actual
print(f"\nTotal Load Power = {P_load1_actual:.4f} + {P_load2_actual:.4f} + {P_load3_actual:.4f} + {P_load4_actual:.4f} + {P_load5_actual:.4f} = {total_load_power:.4f} MW")

print("\n" + "="*50 + "\n")

# Step 2: Calculate total real power generation supplied to the bus
print("Step 2: Calculating Total Real Power Generation\n")

# Hydro Generator
P_hydro_gen = S_hydro_gen * pf_sources
print(f"Hydro Gen. Real Power = {S_hydro_gen:.2f} MVA * {pf_sources:.2f} PF = {P_hydro_gen:.4f} MW")

# Mini Hydro Generator
P_mini_hydro_gen = S_mini_hydro_gen * pf_sources
print(f"Mini Hydro Gen. Real Power = {S_mini_hydro_gen:.2f} MVA * {pf_sources:.2f} PF = {P_mini_hydro_gen:.4f} MW")

# Photovoltaic System (with transmission loss)
P_pv_generated = S_pv * pf_pv
transmission_loss_factor_pv = 1 - pv_transmission_loss_pct / 100
P_pv_delivered = P_pv_generated * transmission_loss_factor_pv
print(f"Photovoltaic Real Power Generated = {S_pv:.2f} MVA * {pf_pv:.2f} PF = {P_pv_generated:.4f} MW")
print(f"Photovoltaic Real Power Delivered (with {pv_transmission_loss_pct}% transmission loss) = {P_pv_generated:.4f} MW * {transmission_loss_factor_pv:.2f} = {P_pv_delivered:.4f} MW")

# Total Generation Power
total_generation_power = P_hydro_gen + P_mini_hydro_gen + P_pv_delivered
print(f"\nTotal Generation Power = {P_hydro_gen:.4f} + {P_mini_hydro_gen:.4f} + {P_pv_delivered:.4f} = {total_generation_power:.4f} MW")

print("\n" + "="*50 + "\n")

# Step 3: Calculate the net real power demand on Bus 11
print("Step 3: Calculating Net Real Power Demand on Bus 11\n")
net_demand = total_load_power - total_generation_power

print("Net Demand = Total Load Power - Total Generation Power")
print(f"Net Demand = ({P_load1_actual:.4f} + {P_load2_actual:.4f} + {P_load3_actual:.4f} + {P_load4_actual:.4f} + {P_load5_actual:.4f}) - ({P_hydro_gen:.4f} + {P_mini_hydro_gen:.4f} + {P_pv_delivered:.4f})")
print(f"Net Demand = {total_load_power:.4f} MW - {total_generation_power:.4f} MW")
print(f"Final Net Real Power Demand = {net_demand:.4f} MW")

print("\n<<<4.2929>>>")