import math

# Step 1: Define initial values from the problem description
# Load values (MW)
p_load1 = 1.89
p_load2_base = 1.75
p_load3 = 1.7
p_load4_base = 1.9
p_load5 = 2.4

# Generation values
s_hydro = 2.0  # MVA
s_mini_hydro = 1.2  # MVA
s_pv = 3.0  # MVA

# Power factors
pf_sources = 0.9
pf_pv = 0.92

# Losses
harmonic_loss_factor = 1.05  # 5% loss means total is 105%
pv_transmission_loss_factor = 1 - 0.04 # 4% loss

# Step 2: Calculate actual load powers including losses
p_load2_actual = p_load2_base * harmonic_loss_factor
p_load4_actual = p_load4_base * harmonic_loss_factor

# Step 3: Calculate total real power load
total_load = p_load1 + p_load2_actual + p_load3 + p_load4_actual + p_load5

# Step 4: Calculate real power generation delivered to the bus
p_hydro_gen = s_hydro * pf_sources
p_mini_hydro_gen = s_mini_hydro * pf_sources
p_pv_gen_gross = s_pv * pf_pv
p_pv_gen_delivered = p_pv_gen_gross * pv_transmission_loss_factor

# Step 5: Calculate total real power generation
total_generation = p_hydro_gen + p_mini_hydro_gen + p_pv_gen_delivered

# Step 6: Calculate the net real power demand on Bus 11
net_demand = total_load - total_generation

# Step 7: Print the detailed equation and the final result
print("Calculation of Net Real Power Demand on Bus 11 (in MW):")
print("\nStep 1: Calculate Total Load")
print(f"Total Load = Load 1 + Load 2 (incl. loss) + Load 3 + Load 4 (incl. loss) + Load 5")
print(f"Total Load = {p_load1} + ({p_load2_base} * {harmonic_loss_factor}) + {p_load3} + ({p_load4_base} * {harmonic_loss_factor}) + {p_load5}")
print(f"Total Load = {p_load1:.2f} + {p_load2_actual:.4f} + {p_load3:.1f} + {p_load4_actual:.3f} + {p_load5:.1f} = {total_load:.4f} MW")

print("\nStep 2: Calculate Total Generation")
print(f"Total Generation = Hydro Gen + Mini Hydro Gen + PV Gen (net of loss)")
print(f"Total Generation = ({s_hydro} * {pf_sources}) + ({s_mini_hydro} * {pf_sources}) + (({s_pv} * {pf_pv}) * {pv_transmission_loss_factor})")
print(f"Total Generation = {p_hydro_gen:.1f} + {p_mini_hydro_gen:.2f} + {p_pv_gen_delivered:.4f} = {total_generation:.4f} MW")

print("\nStep 3: Calculate Net Demand")
print(f"Net Demand = Total Load - Total Generation")
print(f"Net Demand = {total_load:.4f} - {total_generation:.4f}")
print(f"Net Demand = {net_demand:.4f} MW")

# Final equation with all numbers
print("\nFinal Equation:")
print(f"Net Demand = ({p_load1} + {p_load2_actual:.4f} + {p_load3} + {p_load4_actual:.3f} + {p_load5}) - ({p_hydro_gen:.1f} + {p_mini_hydro_gen:.2f} + {p_pv_gen_delivered:.4f})")
print(f"Net Demand = {total_load:.4f} - {total_generation:.4f} = {net_demand:.4f} MW")

# Final answer in the required format
# print(f'<<<{net_demand:.4f}>>>')