import math

def calculate_net_power_demand():
    """
    Calculates the total net real power demand on Bus 11 based on the microgrid diagram and given conditions.
    """
    # --- Step 1: Calculate Total Load Power ---

    # Base real power values for loads from the diagram (in MW)
    p_load1 = 1.89
    p_load2 = 1.75
    p_load3 = 1.7
    p_load4 = 1.9
    p_load5 = 2.4

    # Apply 5% harmonic power loss to Load 2 and Load 4
    harmonic_loss_multiplier = 1 + 0.05
    p_load2_effective = p_load2 * harmonic_loss_multiplier
    p_load4_effective = p_load4 * harmonic_loss_multiplier
    
    # Sum of all effective load powers
    total_load_power = p_load1 + p_load2_effective + p_load3 + p_load4_effective + p_load5

    # --- Step 2: Calculate Total Generation Power ---

    # Apparent power (MVA) and power factors for sources
    s_hydro_gen = 2.0         # MVA
    pf_hydro_gen = 0.9        # All sources except PV have pf=0.9
    
    s_mini_hydro_gen = 1.2    # MVA
    pf_mini_hydro_gen = 0.9
    
    s_pv_gen = 3.0            # MVA
    pf_pv = 0.92

    # Calculate real power (MW) from apparent power (MVA) and power factor (PF)
    # P = S * PF
    p_hydro_gen = s_hydro_gen * pf_hydro_gen
    p_mini_hydro_gen = s_mini_hydro_gen * pf_mini_hydro_gen
    
    # Calculate PV power generated before loss
    p_pv_generated = s_pv_gen * pf_pv
    
    # Apply 4% transmission loss to PV system
    transmission_loss_factor = 1 - 0.04
    p_pv_delivered = p_pv_generated * transmission_loss_factor
    
    # Sum of all delivered generation powers
    total_generation_power = p_hydro_gen + p_mini_hydro_gen + p_pv_delivered

    # --- Step 3: Calculate Net Power Demand ---
    
    net_demand = total_load_power - total_generation_power

    # --- Step 4: Output the results ---
    
    print("Calculating the net real power demand on Bus 11 (all values in MW)\n")

    print("Total Load Power = P_L1 + P_L2_eff + P_L3 + P_L4_eff + P_L5")
    print(f"Total Load Power = {p_load1:.2f} + ({p_load2:.2f} * {harmonic_loss_multiplier:.2f}) + {p_load3:.1f} + ({p_load4:.2f} * {harmonic_loss_multiplier:.2f}) + {p_load5:.1f}")
    print(f"Total Load Power = {p_load1:.2f} + {p_load2_effective:.4f} + {p_load3:.1f} + {p_load4_effective:.4f} + {p_load5:.1f} = {total_load_power:.4f} MW\n")
    
    print("Total Generation Power = P_Hydro + P_MiniHydro + P_PV_del")
    print(f"Total Generation Power = ({s_hydro_gen:.1f} * {pf_hydro_gen:.1f}) + ({s_mini_hydro_gen:.1f} * {pf_mini_hydro_gen:.1f}) + (({s_pv_gen:.1f} * {pf_pv:.2f}) * {transmission_loss_factor:.2f})")
    print(f"Total Generation Power = {p_hydro_gen:.2f} + {p_mini_hydro_gen:.2f} + {p_pv_delivered:.4f} = {total_generation_power:.4f} MW\n")

    print("Final Equation for Net Demand:")
    print(f"Net Demand = ({p_load1:.2f} + {p_load2_effective:.4f} + {p_load3:.1f} + {p_load4_effective:.4f} + {p_load5:.1f}) - ({p_hydro_gen:.2f} + {p_mini_hydro_gen:.2f} + {p_pv_delivered:.4f})")
    print(f"Net Demand = {total_load_power:.4f} MW - {total_generation_power:.4f} MW")
    print(f"Total Net Real Power Demand on Bus 11 = {net_demand:.4f} MW")

if __name__ == '__main__':
    calculate_net_power_demand()
    # To directly return the final value for the system
    # net_demand = (1.89 + 1.75*1.05 + 1.7 + 1.9*1.05 + 2.4) - (2.0*0.9 + 1.2*0.9 + 3.0*0.92*0.96)
    # print(f"<<<{net_demand:.4f}>>>")