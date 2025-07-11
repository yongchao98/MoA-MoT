import math

def calculate_net_power_demand():
    """
    Calculates the total net real power demand on Bus 11 based on the provided microgrid data.
    """

    # --- 1. Define constants and given values ---
    
    # Loads (Real Power in MW from diagram)
    p_load1 = 1.89  # MW
    p_load2_nominal = 1.75  # MW
    p_load3 = 1.70  # MW
    p_load4_nominal = 1.90  # MW
    p_load5 = 2.40  # MW

    # Losses
    harmonic_loss_factor = 1.05  # 5% harmonic loss -> 1 + 0.05
    pv_transmission_efficiency = 0.96 # 4% transmission loss -> 1 - 0.04
    
    # Generators (Apparent Power in MVA from diagram)
    s_hydro = 2.0    # MVA
    s_mini_hydro = 1.2 # MVA
    s_pv = 3.0       # MVA

    # Power Factors
    pf_sources = 0.90  # For Hydro and Mini Hydro
    pf_pv = 0.92       # For PV

    # --- 2. Calculate Total Load Power ---
    
    # Adjust loads with harmonic loss
    p_load2_actual = p_load2_nominal * harmonic_loss_factor
    p_load4_actual = p_load4_nominal * harmonic_loss_factor

    # Sum of all loads
    total_load_power = p_load1 + p_load2_actual + p_load3 + p_load4_actual + p_load5

    # --- 3. Calculate Total Generation Power ---
    
    # Calculate real power from each generator
    p_hydro = s_hydro * pf_sources
    p_mini_hydro = s_mini_hydro * pf_sources
    
    # Calculate PV power delivered to the bus (after loss)
    p_pv_delivered = s_pv * pf_pv * pv_transmission_efficiency

    # Sum of all generation
    total_generation_power = p_hydro + p_mini_hydro + p_pv_delivered

    # --- 4. Calculate Net Real Power Demand ---
    
    net_demand = total_load_power - total_generation_power

    # --- 5. Print the output as a detailed equation ---

    print("The net real power demand on Bus 11 is calculated as: (Total Loads) - (Total Generation)\n")

    # Final Equation with all numbers
    equation_str = (
        f"Net Demand (MW) = \n"
        f"  (Load1 + Load2_actual + Load3 + Load4_actual + Load5) \n"
        f"- (P_Hydro + P_MiniHydro + P_PV_delivered)\n\n"
        f"Net Demand (MW) = \n"
        f"  ({p_load1:.2f} + ({p_load2_nominal:.2f} * {harmonic_loss_factor:.2f}) + {p_load3:.2f} + ({p_load4_nominal:.2f} * {harmonic_loss_factor:.2f}) + {p_load5:.2f}) \n"
        f"- (({s_hydro:.1f} * {pf_sources:.2f}) + ({s_mini_hydro:.1f} * {pf_sources:.2f}) + ({s_pv:.1f} * {pf_pv:.2f} * {pv_transmission_efficiency:.2f}))\n\n"
        f"Net Demand (MW) = \n"
        f"  ({p_load1:.2f} + {p_load2_actual:.4f} + {p_load3:.2f} + {p_load4_actual:.4f} + {p_load5:.2f}) \n"
        f"- ({p_hydro:.2f} + {p_mini_hydro:.2f} + {p_pv_delivered:.4f})\n\n"
        f"Net Demand (MW) = {total_load_power:.4f} - {total_generation_power:.4f}\n"
        f"Net Demand (MW) = {net_demand:.4f}"
    )
    
    print(equation_str)
    
    # Return the final numerical value for the grading format
    return net_demand

# Execute the calculation and print the final result in the required format
final_answer = calculate_net_power_demand()
print(f"\n<<< {final_answer:.4f} >>>")
