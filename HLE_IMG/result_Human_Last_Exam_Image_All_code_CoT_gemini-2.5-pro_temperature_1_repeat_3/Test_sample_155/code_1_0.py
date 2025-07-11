def calculate_net_power_demand():
    """
    Calculates the total net real power demand on Bus 11 based on the provided microgrid data.
    """

    # --- Step 1: Define all given parameters ---

    # Load data (Real Power in MW)
    p_load_1 = 1.89
    p_load_2 = 1.75
    p_load_3 = 1.7
    p_load_4 = 1.9
    p_load_5 = 2.4

    # Generator data (Apparent Power in MVA)
    s_hydro_gen = 2.0
    s_mini_hydro_gen = 1.2
    s_pv_gen = 3.0

    # Power factors
    pf_sources = 0.9  # for all sources except PV
    pf_pv = 0.92      # for PV system

    # Losses
    harmonic_loss_pct = 5.0
    pv_trans_loss_pct = 4.0

    # --- Step 2: Calculate total real power load ---
    p_load_2_eff = p_load_2 * (1 + harmonic_loss_pct / 100)
    p_load_4_eff = p_load_4 * (1 + harmonic_loss_pct / 100)

    total_load_power = p_load_1 + p_load_2_eff + p_load_3 + p_load_4_eff + p_load_5

    print("--- Calculating Total Load Power ---")
    print(f"Total Load Power = P_L1 + P_L2_eff + P_L3 + P_L4_eff + P_L5")
    print(f"Total Load Power = {p_load_1:.2f} MW + ({p_load_2:.2f} * 1.05) MW + {p_load_3:.2f} MW + ({p_load_4:.2f} * 1.05) MW + {p_load_5:.2f} MW")
    print(f"Total Load Power = {p_load_1:.4f} + {p_load_2_eff:.4f} + {p_load_3:.4f} + {p_load_4_eff:.4f} + {p_load_5:.4f} = {total_load_power:.4f} MW\n")

    # --- Step 3: Calculate total real power generation ---
    p_hydro_gen = s_hydro_gen * pf_sources
    p_mini_hydro_gen = s_mini_hydro_gen * pf_sources
    p_pv_generated = s_pv_gen * pf_pv
    p_pv_injected = p_pv_generated * (1 - pv_trans_loss_pct / 100)

    total_generation_power = p_hydro_gen + p_mini_hydro_gen + p_pv_injected

    print("--- Calculating Total Generation Power ---")
    print(f"Total Generation = P_Hydro + P_MiniHydro + P_PV_net")
    print(f"Total Generation = ({s_hydro_gen:.2f} MVA * {pf_sources}) + ({s_mini_hydro_gen:.2f} MVA * {pf_sources}) + (({s_pv_gen:.2f} MVA * {pf_pv}) * (1 - {pv_trans_loss_pct/100}))")
    print(f"Total Generation Power = {p_hydro_gen:.4f} + {p_mini_hydro_gen:.4f} + {p_pv_injected:.4f} = {total_generation_power:.4f} MW\n")

    # --- Step 4: Calculate net real power demand ---
    net_demand = total_load_power - total_generation_power

    print("--- Calculating Net Real Power Demand on Bus 11 ---")
    print(f"Net Demand = Total Load Power - Total Generation Power")
    print(f"Net Demand = {total_load_power:.4f} MW - {total_generation_power:.4f} MW = {net_demand:.4f} MW")

if __name__ == '__main__':
    calculate_net_power_demand()