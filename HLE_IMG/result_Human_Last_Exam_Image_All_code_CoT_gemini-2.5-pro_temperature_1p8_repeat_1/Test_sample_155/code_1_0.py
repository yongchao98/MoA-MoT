import math

# --- Step 1: Calculate the total real power load on Bus 11 ---
print("--- Step 1: Calculating Total Real Power Load ---")

# Load 1
load1_p = 1.89  # MW
print(f"Load 1 Real Power = {load1_p:.2f} MW")

# Load 2
load2_p_base = 1.75  # MW
harmonic_loss_pct = 5
load2_p_adj = load2_p_base * (1 + harmonic_loss_pct / 100)
print(f"Load 2 Real Power = {load2_p_base:.2f} MW, with a {harmonic_loss_pct}% harmonic loss. Adjusted Power = {load2_p_adj:.4f} MW")

# Load 3
load3_p = 1.7  # MW
print(f"Load 3 Real Power = {load3_p:.2f} MW (Power factor is not needed for real power)")

# Load 4
load4_p_base = 1.9  # MW
load4_p_adj = load4_p_base * (1 + harmonic_loss_pct / 100)
print(f"Load 4 Real Power = {load4_p_base:.2f} MW, with a {harmonic_loss_pct}% harmonic loss. Adjusted Power = {load4_p_adj:.4f} MW")

# Load 5
load5_p = 2.4 # MW
print(f"Load 5 Real Power = {load5_p:.2f} MW (Additional reactive power does not affect real power)")

# Calculate Total Load
total_load_p = load1_p + load2_p_adj + load3_p + load4_p_adj + load5_p
print("\nCalculating Total Load:")
print(f"Total Load = {load1_p:.2f} + {load2_p_adj:.4f} + {load3_p:.2f} + {load4_p_adj:.4f} + {load5_p:.2f} = {total_load_p:.4f} MW")

print("\n" + "="*50 + "\n")

# --- Step 2: Calculate the total real power generation supplied to Bus 11 ---
print("--- Step 2: Calculating Total Real Power Generation ---")

# General source Power Factor
source_pf = 0.9

# Hydro Generator
hydro_gen_s = 2.0  # MVA
hydro_gen_p = hydro_gen_s * source_pf
print(f"Hydro Gen. Real Power = {hydro_gen_s:.1f} MVA * {source_pf} PF = {hydro_gen_p:.2f} MW")

# Mini Hydro Generator
mini_hydro_gen_s = 1.2 # MVA
mini_hydro_gen_p = mini_hydro_gen_s * source_pf
print(f"Mini Hydro Gen. Real Power = {mini_hydro_gen_s:.1f} MVA * {source_pf} PF = {mini_hydro_gen_p:.2f} MW")

# Photovoltaic System
pv_s = 3.0  # MVA
pv_pf = 0.92
pv_transmission_loss_pct = 4
pv_p_initial = pv_s * pv_pf
pv_p_delivered = pv_p_initial * (1 - pv_transmission_loss_pct / 100)
print(f"PV Gen. Initial Real Power = {pv_s:.1f} MVA * {pv_pf} PF = {pv_p_initial:.2f} MW")
print(f"PV Gen. Delivered Power after {pv_transmission_loss_pct}% loss = {pv_p_initial:.2f} MW * (1 - {pv_transmission_loss_pct/100}) = {pv_p_delivered:.4f} MW")

# Calculate Total Generation
total_gen_p = hydro_gen_p + mini_hydro_gen_p + pv_p_delivered
print("\nCalculating Total Generation:")
print(f"Total Generation = {hydro_gen_p:.2f} + {mini_hydro_gen_p:.2f} + {pv_p_delivered:.4f} = {total_gen_p:.4f} MW")

print("\n" + "="*50 + "\n")

# --- Step 3: Calculate the total net real power demand on Bus 11 ---
print("--- Step 3: Calculating Net Real Power Demand on Bus 11 ---")
net_demand_p = total_load_p - total_gen_p

print("Net Demand = Total Load - Total Generation")

print("\nFinal Equation using all numbers:")
print(f"Net Demand (MW) = ({load1_p} + {load2_p_base} * (1 + {harmonic_loss_pct/100}) + {load3_p} + {load4_p_base} * (1 + {harmonic_loss_pct/100}) + {load5_p}) - (({hydro_gen_s} * {source_pf}) + ({mini_hydro_gen_s} * {source_pf}) + ({pv_s} * {pv_pf} * (1 - {pv_transmission_loss_pct/100})))")

print(f"\nNet Demand = {total_load_p:.4f} MW - {total_gen_p:.4f} MW")
print(f"Net Demand = {net_demand_p:.4f} MW")

print(f"\n<<<{net_demand_p:.4f}>>>")