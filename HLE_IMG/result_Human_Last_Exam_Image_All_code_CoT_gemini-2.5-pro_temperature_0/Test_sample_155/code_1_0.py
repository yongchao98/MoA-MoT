import math

# --- Input Values from the Problem ---

# Loads (in MW)
load1_p = 1.89
load2_p = 1.75
load3_p = 1.7
load4_p = 1.9
load5_p = 2.4

# Generators (Apparent Power in MVA)
hydro_gen_s = 2.0
mini_hydro_gen_s = 1.2
pv_gen_s = 3.0

# Power Factors
source_pf = 0.90
pv_pf = 0.92

# Losses
harmonic_loss_pct = 0.05
pv_transmission_loss_pct = 0.04

# --- Calculation Steps ---

# 1. Calculate Total Load Demand
print("Step 1: Calculate Total Load Demand")
# Adjust loads for harmonic loss
load2_p_adj = load2_p * (1 + harmonic_loss_pct)
load4_p_adj = load4_p * (1 + harmonic_loss_pct)
print(f"Adjusted Load 2 (with {harmonic_loss_pct*100}% harmonic loss): {load2_p} MW * {1 + harmonic_loss_pct} = {load2_p_adj:.4f} MW")
print(f"Adjusted Load 4 (with {harmonic_loss_pct*100}% harmonic loss): {load4_p} MW * {1 + harmonic_loss_pct} = {load4_p_adj:.4f} MW")

# Sum all loads
total_load_p = load1_p + load2_p_adj + load3_p + load4_p_adj + load5_p
print(f"Total Load = {load1_p} + {load2_p_adj:.4f} + {load3_p} + {load4_p_adj:.4f} + {load5_p} = {total_load_p:.4f} MW\n")


# 2. Calculate Total Power Generation
print("Step 2: Calculate Total Delivered Power Generation")
# Calculate real power for each generator
hydro_gen_p = hydro_gen_s * source_pf
mini_hydro_gen_p = mini_hydro_gen_s * source_pf
pv_gen_p_initial = pv_gen_s * pv_pf
print(f"Hydro Gen. Real Power = {hydro_gen_s} MVA * {source_pf} = {hydro_gen_p:.4f} MW")
print(f"Mini Hydro Gen. Real Power = {mini_hydro_gen_s} MVA * {source_pf} = {mini_hydro_gen_p:.4f} MW")
print(f"PV System Generated Power = {pv_gen_s} MVA * {pv_pf} = {pv_gen_p_initial:.4f} MW")

# Adjust PV generation for transmission loss
pv_gen_p_delivered = pv_gen_p_initial * (1 - pv_transmission_loss_pct)
print(f"PV Delivered Power (with {pv_transmission_loss_pct*100}% loss) = {pv_gen_p_initial:.4f} MW * (1 - {pv_transmission_loss_pct}) = {pv_gen_p_delivered:.4f} MW")

# Sum all delivered generation
total_generation_p = hydro_gen_p + mini_hydro_gen_p + pv_gen_p_delivered
print(f"Total Generation = {hydro_gen_p:.4f} + {mini_hydro_gen_p:.4f} + {pv_gen_p_delivered:.4f} = {total_generation_p:.4f} MW\n")


# 3. Calculate Net Power Demand
print("Step 3: Calculate Net Real Power Demand on Bus 11")
net_demand_p = total_load_p - total_generation_p
print("Net Demand = Total Load - Total Generation")
print(f"Net Demand = {total_load_p:.4f} MW - {total_generation_p:.4f} MW")
print(f"The total net real power demand on Bus 11 is: {net_demand_p:.4f} MW")

print("\n--- Final Equation ---")
print(f"Net Demand = ({load1_p} + {load2_p}*(1+{harmonic_loss_pct}) + {load3_p} + {load4_p}*(1+{harmonic_loss_pct}) + {load5_p}) - ({hydro_gen_s}*{source_pf} + {mini_hydro_gen_s}*{source_pf} + {pv_gen_s}*{pv_pf}*(1-{pv_transmission_loss_pct}))")
print(f"Net Demand = {total_load_p:.4f} - {total_generation_p:.4f} = {net_demand_p:.4f} MW")
<<<4.2929>>>