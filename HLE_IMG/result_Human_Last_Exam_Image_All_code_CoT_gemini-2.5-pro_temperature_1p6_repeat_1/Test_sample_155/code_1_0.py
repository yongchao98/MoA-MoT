import math

# Step 1: Define base load values from the diagram
load1_p = 1.89  # MW
load2_p = 1.75  # MW
load3_p = 1.7   # MW
load4_p = 1.9   # MW
load5_p = 2.4   # MW

# Step 2: Define adjustments and source parameters from the problem description
harmonic_loss_pct = 5.0
pv_s = 3.0      # MVA
pv_pf = 0.92
pv_transmission_loss_pct = 4.0

# Step 3: Calculate adjusted load values
harmonic_loss_multiplier = 1 + (harmonic_loss_pct / 100)
load2_p_adj = load2_p * harmonic_loss_multiplier
load4_p_adj = load4_p * harmonic_loss_multiplier

# Step 4: Calculate total real power load
total_load_p = load1_p + load2_p_adj + load3_p + load4_p_adj + load5_p

# Step 5: Calculate delivered real power from the PV system
pv_p_generated = pv_s * pv_pf
pv_transmission_loss = pv_p_generated * (pv_transmission_loss_pct / 100)
pv_p_delivered = pv_p_generated - pv_transmission_loss

# Step 6: Calculate the net real power demand on Bus 11
net_demand_p = total_load_p - pv_p_delivered

# Step 7: Print the calculation steps
print("Calculating Total Real Power Load:")
print(f"Load 1 = {load1_p:.2f} MW")
print(f"Load 2 (adjusted for {harmonic_loss_pct}% harmonic loss) = {load2_p:.2f} MW * {harmonic_loss_multiplier:.2f} = {load2_p_adj:.4f} MW")
print(f"Load 3 = {load3_p:.2f} MW")
print(f"Load 4 (adjusted for {harmonic_loss_pct}% harmonic loss) = {load4_p:.2f} MW * {harmonic_loss_multiplier:.2f} = {load4_p_adj:.4f} MW")
print(f"Load 5 = {load5_p:.2f} MW")
print("-" * 30)
print(f"Total Load = {load1_p:.2f} + {load2_p_adj:.4f} + {load3_p:.2f} + {load4_p_adj:.4f} + {load5_p:.2f} = {total_load_p:.4f} MW")
print("\nCalculating Delivered Real Power from PV System:")
print(f"PV Power Generated = {pv_s:.2f} MVA * {pv_pf:.2f} PF = {pv_p_generated:.4f} MW")
print(f"PV Power Delivered to Bus (after {pv_transmission_loss_pct}% loss) = {pv_p_generated:.4f} MW * (1 - {pv_transmission_loss_pct/100:.2f}) = {pv_p_delivered:.4f} MW")
print("\nCalculating Net Real Power Demand on Bus 11:")
print(f"Net Demand = Total Load - PV Power Delivered")
print(f"Net Demand = {total_load_p:.4f} MW - {pv_p_delivered:.4f} MW = {net_demand_p:.4f} MW")

# Final answer
print(f"\nThe total net real power demand on Bus 11 is {net_demand_p:.4f} MW.")
print(f"<<<{net_demand_p:.4f}>>>")