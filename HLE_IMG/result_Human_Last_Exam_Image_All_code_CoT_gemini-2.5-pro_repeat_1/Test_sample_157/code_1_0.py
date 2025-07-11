import cmath
import math

def calculate_power_losses():
    """
    Calculates the total resistive power losses in the power system.
    """
    # 1. Establish Base Values
    S_base_MVA = 100.0  # MVA
    V_base_kV = 230.0   # kV
    Z_base = (V_base_kV**2) / S_base_MVA  # Ohms

    # Given data
    # Transmission line resistances in Ohms
    R_actual = {
        'T7_8': 4.50, 'T8_9': 6.30,
        'T5_7': 16.93, 'T6_9': 20.63,
        'T4_5': 5.29, 'T4_6': 8.99
    }

    # Load data
    P_base_load_MW = 50.0 / 3.0
    P_add_load_MW = 100.0
    load_pf = 0.9
    
    # 2. Calculate Per-Unit (p.u.) Resistances
    R_pu = {line: res / Z_base for line, res in R_actual.items()}

    # 3. Determine Power Injections and Withdrawals in p.u.
    # Calculate reactive power component from power factor
    tan_theta = math.tan(math.acos(load_pf))

    # Calculate load powers in MVA, then convert to p.u.
    S_L5_pu = (P_base_load_MW + 1j * P_base_load_MW * tan_theta) / S_base_MVA
    S_L6_pu = (P_base_load_MW + 1j * P_base_load_MW * tan_theta) / S_base_MVA
    P_L8_total_MW = P_base_load_MW + P_add_load_MW
    S_L8_pu = (P_L8_total_MW + 1j * P_L8_total_MW * tan_theta) / S_base_MVA

    # Power Generation Dispatch Assumption
    P_total_load_MW = 3 * P_base_load_MW + P_add_load_MW
    P_gen_each_MW = P_total_load_MW / 3.0
    S_G2_pu = P_gen_each_MW / S_base_MVA  # Assume Q=0 for VSG 2
    S_G3_pu = P_gen_each_MW / S_base_MVA  # Assume Q=0 for VSG 3

    # 4. Approximate Power Flow in p.u.
    # Assumption: Load at Bus 8 is split equally between T7_8 and T8_9
    S_flow_78 = S_L8_pu / 2.0
    S_flow_98 = S_L8_pu / 2.0

    # Calculate other flows using Kirchhoff's Current Law for power (S)
    # At Bus 7: S_G2_in = S_78_out + S_75_out => S_75 = S_G2 - S_78
    S_flow_75 = S_G2_pu - S_flow_78
    S_flow_57 = -S_flow_75

    # At Bus 9: S_G3_in = S_98_out + S_96_out => S_96 = S_G3 - S_98
    S_flow_96 = S_G3_pu - S_flow_98
    S_flow_69 = -S_flow_96

    # At Bus 5: S_45_in + S_75_in = S_L5_out => S_45 = S_L5 + S_57
    S_flow_45 = S_L5_pu + S_flow_57

    # At Bus 6: S_46_in + S_96_in = S_L6_out => S_46 = S_L6 + S_69
    S_flow_46 = S_L6_pu + S_flow_69

    # Store flows for calculation. Magnitudes will be used.
    S_flows_pu = {
        'T7_8': S_flow_78, 'T8_9': S_flow_98, # Note: T8_9 flow is S_98
        'T5_7': S_flow_57, 'T6_9': S_flow_69,
        'T4_5': S_flow_45, 'T4_6': S_flow_46
    }

    # 5. Calculate Power Losses
    # P_loss (MW) = |S_pu|^2 * R_pu * S_base_MVA
    losses_MW = {}
    total_loss_MW = 0
    
    print("Calculating losses for each transmission line:")
    
    final_equation_terms = []
    
    for line, S_flow in S_flows_pu.items():
        # Using abs(S_flow)**2 as an approximation for |I_pu|^2
        loss_pu = (abs(S_flow)**2) * R_pu[line]
        loss_MW = loss_pu * S_base_MVA
        losses_MW[line] = loss_MW
        total_loss_MW += loss_MW
        
        print(f"Loss in {line}: {loss_MW:.3f} MW")
        final_equation_terms.append(f"{loss_MW:.3f}")

    print("\nTotal Power Loss Calculation:")
    final_equation = "Total Loss = " + " + ".join(final_equation_terms)
    print(f"{final_equation} = {total_loss_MW:.3f} MW")

    return total_loss_MW

# Execute the calculation and print the final result
total_losses = calculate_power_losses()
print(f"\nThe total power losses in the system are {total_losses:.3f} MW.")
print(f"<<<{total_losses:.3f}>>>")
