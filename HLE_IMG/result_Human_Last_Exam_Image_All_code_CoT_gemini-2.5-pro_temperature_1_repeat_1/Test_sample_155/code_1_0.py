import math

# Step 1: Define all the given values from the problem description and diagram.

# Load values (in MW)
p_load1 = 1.89
p_load2_base = 1.75
p_load3 = 1.7
p_load4_base = 1.9
p_load5 = 2.4

# Generator values
s_hydro_gen = 2.0  # MVA
s_mini_hydro_gen = 1.2  # MVA
s_pv = 3.0  # MVA

# Power factors
pf_sources = 0.9  # For Hydro and Mini Hydro
pf_pv = 0.92      # For Photovoltaic system

# Losses
harmonic_loss_percent = 5.0
pv_transmission_loss_percent = 4.0

# Step 2: Calculate the total load demand, including losses.

# Harmonic loss increases the power drawn by the load
harmonic_loss_multiplier = 1 + (harmonic_loss_percent / 100)
p_load2_actual = p_load2_base * harmonic_loss_multiplier
p_load4_actual = p_load4_base * harmonic_loss_multiplier

# Sum of all loads
total_load_demand = p_load1 + p_load2_actual + p_load3 + p_load4_actual + p_load5

# Step 3: Calculate the total power generation, including losses.

# Real power from Hydro and Mini Hydro generators
p_hydro_gen = s_hydro_gen * pf_sources
p_mini_hydro_gen = s_mini_hydro_gen * pf_sources

# Real power from PV system, accounting for transmission loss
pv_transmission_efficiency = 1 - (pv_transmission_loss_percent / 100)
p_pv_delivered = s_pv * pf_pv * pv_transmission_efficiency

# Sum of all generation
total_generation = p_hydro_gen + p_mini_hydro_gen + p_pv_delivered

# Step 4: Calculate the net real power demand on Bus 11.
net_demand = total_load_demand - total_generation

# Step 5: Print the final equation and the result.
print("The net real power demand on Bus 11 is calculated as: (Total Load) - (Total Generation)")
print("\nFinal Equation (in MW):")

# Constructing and printing the detailed equation string
equation = (
    f"Net Demand = "
    f"({p_load1} + {p_load2_base} * {harmonic_loss_multiplier} + {p_load3} + {p_load4_base} * {harmonic_loss_multiplier} + {p_load5}) - "
    f"({s_hydro_gen} * {pf_sources} + {s_mini_hydro_gen} * {pf_sources} + {s_pv} * {pf_pv} * {pv_transmission_efficiency})"
)
print(equation)

# Printing the calculated values and the final result
calculation_summary = (
    f"Net Demand = "
    f"({p_load1} + {p_load2_actual:.4f} + {p_load3} + {p_load4_actual:.4f} + {p_load5}) - "
    f"({p_hydro_gen:.4f} + {p_mini_hydro_gen:.4f} + {p_pv_delivered:.4f})"
)
print(calculation_summary)

result_summary = (
    f"Net Demand = {total_load_demand:.4f} - {total_generation:.4f} = {net_demand:.4f} MW"
)
print(result_summary)

print(f"\n<<< {net_demand:.4f} >>>")