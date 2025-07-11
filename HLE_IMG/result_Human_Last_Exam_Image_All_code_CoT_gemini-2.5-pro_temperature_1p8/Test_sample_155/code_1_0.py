import math

def calculate_net_power_demand():
    """
    Calculates the total net real power demand on Bus 11 based on the microgrid diagram and given conditions.
    """
    
    # --- Given Data ---

    # Loads (Real Power in MW from the diagram)
    p_load1 = 1.89  # MW
    p_load2 = 1.75  # MW
    p_load3 = 1.7   # MW
    p_load4 = 1.9   # MW
    p_load5 = 2.4   # MW

    # Generators (Apparent Power in MVA and Power Factor from the diagram and text)
    s_hydro_gen = 2.0         # MVA
    s_mini_hydro_gen = 1.2    # MVA
    s_pv = 3.0                # MVA
    pf_sources = 0.9          # for Hydro and Mini-Hydro
    pf_pv = 0.92              # for PV

    # Special Conditions (Losses)
    harmonic_loss_factor = 0.05       # 5% for Load 2 and 4
    pv_transmission_loss_factor = 0.04 # 4% for PV system

    # --- Calculations ---

    # 1. Calculate total load demand, including harmonic losses
    p_load2_adj = p_load2 * (1 + harmonic_loss_factor)
    p_load4_adj = p_load4 * (1 + harmonic_loss_factor)
    
    total_load_power = p_load1 + p_load2_adj + p_load3 + p_load4_adj + p_load5

    # 2. Calculate total generated power delivered to Bus 11
    # PV power with transmission loss
    p_pv_generated = s_pv * pf_pv
    p_pv_delivered = p_pv_generated * (1 - pv_transmission_loss_factor)

    # Hydro and Mini-Hydro power (no losses specified)
    p_hydro_delivered = s_hydro_gen * pf_sources
    p_mini_hydro_delivered = s_mini_hydro_gen * pf_sources
    
    total_generated_power = p_pv_delivered + p_hydro_delivered + p_mini_hydro_delivered

    # 3. Calculate the net real power demand on Bus 11
    net_demand = total_load_power - total_generated_power

    # --- Output Results ---
    
    print("The calculation for the total net real power demand on Bus 11 is as follows:\n")

    print("Step 1: Calculate the total power demand from the loads, including harmonic losses.")
    print("P_loads = P_L1 + P_L2_adjusted + P_L3 + P_L4_adjusted + P_L5")
    print(f"P_loads = {p_load1} MW + ({p_load2} MW * 1.05) + {p_load3} MW + ({p_load4} MW * 1.05) + {p_load5} MW")
    print(f"P_loads = {p_load1} MW + {p_load2_adj:.4f} MW + {p_load3} MW + {p_load4_adj:.4f} MW + {p_load5} MW = {total_load_power:.4f} MW\n")

    print("Step 2: Calculate the total power supplied by the generators to Bus 11.")
    print("P_gens = P_PV_delivered + P_Hydro_delivered + P_MiniHydro_delivered")
    print(f"P_gens = ({s_pv} MVA * {pf_pv} * (1-{pv_transmission_loss_factor})) + ({s_hydro_gen} MVA * {pf_sources}) + ({s_mini_hydro_gen} MVA * {pf_sources})")
    print(f"P_gens = {p_pv_delivered:.4f} MW + {p_hydro_delivered:.4f} MW + {p_mini_hydro_delivered:.4f} MW = {total_generated_power:.4f} MW\n")

    print("Step 3: Calculate the net demand on Bus 11 (Total Loads - Total Generation).")
    print("Net Demand = P_loads - P_gens")
    print(f"Net Demand = {total_load_power:.4f} MW - {total_generated_power:.4f} MW")
    print(f"The total net real power demand on Bus 11 is {net_demand:.4f} MW.")
    
    # Final answer in the required format
    print(f"\n<<<{net_demand:.4f}>>>")

# Run the calculation
calculate_net_power_demand()