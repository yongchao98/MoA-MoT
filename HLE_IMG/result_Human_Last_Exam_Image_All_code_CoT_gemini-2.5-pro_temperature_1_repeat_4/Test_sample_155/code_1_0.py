import math

# --- System Parameters ---

# Loads (Base Real Power in MW)
p_load1 = 1.89
p_load2 = 1.75
p_load3 = 1.7
p_load4 = 1.9
p_load5 = 2.4

# Generators (Apparent Power in MVA)
s_hydro_gen = 2.0
s_mini_hydro_gen = 1.2
s_pv = 3.0

# Power Factors
pf_sources = 0.9  # For Hydro Gen, Mini Hydro Gen
pf_pv = 0.92      # For Photovoltaic system

# Losses
harmonic_loss_factor = 1.05  # 5% loss -> demand increases by 5%
pv_transmission_loss_factor = 1 - 0.04 # 4% loss

# --- Calculations ---

# 1. Calculate actual load demand considering losses
p_load2_actual = p_load2 * harmonic_loss_factor
p_load4_actual = p_load4 * harmonic_loss_factor
total_load_demand = p_load1 + p_load2_actual + p_load3 + p_load4_actual + p_load5

# 2. Calculate actual power generation considering power factors and losses
p_hydro_gen = s_hydro_gen * pf_sources
p_mini_hydro_gen = s_mini_hydro_gen * pf_sources
p_pv_net = s_pv * pf_pv * pv_transmission_loss_factor
total_generation = p_hydro_gen + p_mini_hydro_gen + p_pv_net

# 3. Calculate the net real power demand on Bus 11
net_demand = total_load_demand - total_generation

# --- Output the Results ---
print("Calculating the total net real power demand on Bus 11.")
print("\nStep 1: Calculate Total Power Demand from all loads.")
print(f"Total Demand = Load1 + Load2*(1+HarmonicLoss) + Load3 + Load4*(1+HarmonicLoss) + Load5")
print(f"Total Demand = {p_load1} + ({p_load2} * {harmonic_loss_factor}) + {p_load3} + ({p_load4} * {harmonic_loss_factor}) + {p_load5}")
print(f"Total Demand = {p_load1} + {p_load2_actual:.4f} + {p_load3} + {p_load4_actual:.3f} + {p_load5} = {total_load_demand:.4f} MW")

print("\nStep 2: Calculate Total Power Supplied by all generators.")
print(f"Total Generation = HydroGen + MiniHydroGen + PV*(1-TransmissionLoss)")
print(f"Total Generation = ({s_hydro_gen} * {pf_sources}) + ({s_mini_hydro_gen} * {pf_sources}) + (({s_pv} * {pf_pv}) * {pv_transmission_loss_factor})")
print(f"Total Generation = {p_hydro_gen:.2f} + {p_mini_hydro_gen:.2f} + {p_pv_net:.4f} = {total_generation:.4f} MW")

print("\nStep 3: Calculate Net Demand on Bus 11.")
print(f"Net Demand = Total Demand - Total Generation")
print(f"Net Demand = {total_load_demand:.4f} MW - {total_generation:.4f} MW")
print(f"The total net real power demand on Bus 11 is {net_demand:.4f} MW.")

print("\n--- Final Equation ---")
print(f"({p_load1} + {p_load2_actual:.4f} + {p_load3} + {p_load4_actual:.3f} + {p_load5}) - ({p_hydro_gen:.1f} + {p_mini_hydro_gen:.2f} + {p_pv_net:.4f}) = {net_demand:.4f} MW")
<<<4.2929>>>