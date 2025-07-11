import math

def calculate_net_power_demand():
    """
    Calculates the total net real power demand on Bus 11.
    """
    # --- Input Data ---
    # Loads (MW)
    p_load1 = 1.89
    p_load2_base = 1.75
    p_load3 = 1.7
    p_load4_base = 1.9
    p_load5 = 2.4

    # Generators (MVA)
    s_hydro = 2.0
    s_mini_hydro = 1.2
    s_pv = 3.0

    # Power Factors
    pf_sources = 0.9  # for hydro, mini hydro, and load 3
    pf_pv = 0.92

    # Losses and other conditions
    harmonic_loss_pct = 5.0
    pv_transmission_loss_pct = 4.0

    # --- Calculations ---

    # 1. Calculate real power from generators
    print("--- Real Power Generation Calculation ---")
    # Hydro Gen
    p_hydro = s_hydro * pf_sources
    print(f"Hydro Generator Power = {s_hydro:.2f} MVA * {pf_sources:.2f} PF = {p_hydro:.4f} MW")

    # Mini Hydro Gen
    p_mini_hydro = s_mini_hydro * pf_sources
    print(f"Mini Hydro Generator Power = {s_mini_hydro:.2f} MVA * {pf_sources:.2f} PF = {p_mini_hydro:.4f} MW")

    # Photovoltaic Gen
    p_pv_initial = s_pv * pf_pv
    p_pv_delivered = p_pv_initial * (1 - pv_transmission_loss_pct / 100.0)
    print(f"Photovoltaic Power (generated) = {s_pv:.2f} MVA * {pf_pv:.2f} PF = {p_pv_initial:.4f} MW")
    print(f"Photovoltaic Power (delivered to Bus 11 after {pv_transmission_loss_pct}% loss) = {p_pv_initial:.4f} MW * (1 - {pv_transmission_loss_pct/100:.2f}) = {p_pv_delivered:.4f} MW")

    # Total Generation
    p_gen_total = p_hydro + p_mini_hydro + p_pv_delivered
    print("--------------------------------------------------")
    print(f"Total Real Power Generation = {p_hydro:.4f} MW + {p_mini_hydro:.4f} MW + {p_pv_delivered:.4f} MW = {p_gen_total:.4f} MW")
    print("\n")


    # 2. Calculate real power for loads
    print("--- Real Power Load Calculation ---")
    # Load 1
    print(f"Load 1 Power = {p_load1:.4f} MW")
    
    # Load 2 with harmonic loss
    p_load2_actual = p_load2_base * (1 + harmonic_loss_pct / 100.0)
    print(f"Load 2 Power (with {harmonic_loss_pct}% harmonic loss) = {p_load2_base:.2f} MW * (1 + {harmonic_loss_pct/100:.2f}) = {p_load2_actual:.4f} MW")

    # Load 3
    print(f"Load 3 Power = {p_load3:.4f} MW")

    # Load 4 with harmonic loss
    p_load4_actual = p_load4_base * (1 + harmonic_loss_pct / 100.0)
    print(f"Load 4 Power (with {harmonic_loss_pct}% harmonic loss) = {p_load4_base:.2f} MW * (1 + {harmonic_loss_pct/100:.2f}) = {p_load4_actual:.4f} MW")

    # Load 5
    print(f"Load 5 Power = {p_load5:.4f} MW")

    # Total Load
    p_load_total = p_load1 + p_load2_actual + p_load3 + p_load4_actual + p_load5
    print("--------------------------------------------------")
    print(f"Total Real Power Load = {p_load1:.4f} MW + {p_load2_actual:.4f} MW + {p_load3:.4f} MW + {p_load4_actual:.4f} MW + {p_load5:.4f} MW = {p_load_total:.4f} MW")
    print("\n")

    # 3. Calculate Net Real Power Demand on Bus 11
    net_demand = p_load_total - p_gen_total
    print("--- Net Real Power Demand Calculation ---")
    print("Net Demand on Bus 11 = Total Load - Total Generation")
    print(f"Net Demand on Bus 11 = {p_load_total:.4f} MW - {p_gen_total:.4f} MW = {net_demand:.4f} MW")
    
    return net_demand

if __name__ == '__main__':
    final_answer = calculate_net_power_demand()
    print(f"\n<<<The final answer is {final_answer:.4f}>>>")
