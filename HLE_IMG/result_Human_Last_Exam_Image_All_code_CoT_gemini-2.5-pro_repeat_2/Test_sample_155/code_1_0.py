import math

def calculate_net_power_demand():
    """
    Calculates the total net real power demand on Bus 11 based on the provided microgrid diagram and parameters.
    """
    # --- Input Parameters ---

    # Loads (Real Power in MW from the diagram)
    p_load1 = 1.89
    p_load2 = 1.75
    p_load3 = 1.7
    p_load4 = 1.9
    p_load5 = 2.4

    # Generations (Apparent Power in MVA from the diagram)
    s_hydro_gen = 2.0
    s_mini_hydro_gen = 1.2
    s_pv = 3.0

    # Power Factors (PF)
    pf_sources = 0.9      # For Hydro Gen, Mini Hydro Gen, and Load 3
    pf_pv = 0.92          # For PV system

    # Losses
    harmonic_loss_percent = 5.0
    pv_transmission_loss_percent = 4.0

    # --- Calculations ---

    print("Step 1: Calculate Total Real Power Demand from Loads\n")

    # Adjust for harmonic losses
    harmonic_loss_multiplier = 1 + (harmonic_loss_percent / 100)
    p_load2_adj = p_load2 * harmonic_loss_multiplier
    p_load4_adj = p_load4 * harmonic_loss_multiplier

    print(f"Load 2 with {harmonic_loss_percent}% harmonic loss: {p_load2:.2f} MW * {harmonic_loss_multiplier:.2f} = {p_load2_adj:.4f} MW")
    print(f"Load 4 with {harmonic_loss_percent}% harmonic loss: {p_load4:.2f} MW * {harmonic_loss_multiplier:.2f} = {p_load4_adj:.4f} MW\n")

    # Sum of all loads
    total_load_power = p_load1 + p_load2_adj + p_load3 + p_load4_adj + p_load5
    print("Total Load Power = Load 1 + Load 2 (adj) + Load 3 + Load 4 (adj) + Load 5")
    print(f"Total Load Power = {p_load1:.2f} + {p_load2_adj:.4f} + {p_load3:.1f} + {p_load4_adj:.4f} + {p_load5:.1f}")
    print(f"Total Load Power = {total_load_power:.4f} MW\n")
    
    print("--------------------------------------------------\n")

    print("Step 2: Calculate Total Real Power Supplied by Generators\n")

    # Calculate real power from each generator
    p_hydro_gen = s_hydro_gen * pf_sources
    print(f"Hydro Gen real power: {s_hydro_gen:.1f} MVA * {pf_sources:.1f} PF = {p_hydro_gen:.2f} MW")

    p_mini_hydro_gen = s_mini_hydro_gen * pf_sources
    print(f"Mini Hydro Gen real power: {s_mini_hydro_gen:.1f} MVA * {pf_sources:.1f} PF = {p_mini_hydro_gen:.2f} MW")
    
    # Calculate PV power after transmission loss
    p_pv_generated = s_pv * pf_pv
    pv_transmission_loss_multiplier = 1 - (pv_transmission_loss_percent / 100)
    p_pv_supplied = p_pv_generated * pv_transmission_loss_multiplier
    print(f"PV generated real power: {s_pv:.1f} MVA * {pf_pv:.2f} PF = {p_pv_generated:.2f} MW")
    print(f"PV power supplied (after {pv_transmission_loss_percent}% loss): {p_pv_generated:.2f} MW * {pv_transmission_loss_multiplier:.2f} = {p_pv_supplied:.4f} MW\n")

    # Sum of all generation
    total_generation_power = p_hydro_gen + p_mini_hydro_gen + p_pv_supplied
    print("Total Generation Power = Hydro Gen + Mini Hydro Gen + PV (supplied)")
    print(f"Total Generation Power = {p_hydro_gen:.2f} + {p_mini_hydro_gen:.2f} + {p_pv_supplied:.4f}")
    print(f"Total Generation Power = {total_generation_power:.4f} MW\n")

    print("--------------------------------------------------\n")

    print("Step 3: Calculate Net Real Power Demand on Bus 11\n")
    net_demand = total_load_power - total_generation_power
    print("Final Equation: Net Demand = Total Load Power - Total Generation Power")
    print(f"Net Demand = ({p_load1} + {p_load2}*{harmonic_loss_multiplier} + {p_load3} + {p_load4}*{harmonic_loss_multiplier} + {p_load5}) - ({s_hydro_gen}*{pf_sources} + {s_mini_hydro_gen}*{pf_sources} + {s_pv}*{pf_pv}*{pv_transmission_loss_multiplier})")
    print(f"Net Demand = ({p_load1:.2f} + {p_load2_adj:.4f} + {p_load3:.1f} + {p_load4_adj:.4f} + {p_load5:.1f}) - ({p_hydro_gen:.2f} + {p_mini_hydro_gen:.2f} + {p_pv_supplied:.4f})")
    print(f"Net Demand = {total_load_power:.4f} MW - {total_generation_power:.4f} MW")
    print(f"Net Demand = {net_demand:.4f} MW")

    return net_demand

if __name__ == '__main__':
    result = calculate_net_power_demand()
