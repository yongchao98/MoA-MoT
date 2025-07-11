import cmath
import math

def calculate_power_losses():
    """
    Calculates the total resistive power losses in the power system.
    """
    # 1. Define System Parameters
    V_line_kV = 230.0  # Line voltage in kV
    pf = 0.9  # Lagging power factor

    # Base loads (in MW)
    P_load_5 = 50.0 / 3.0
    P_load_6 = 50.0 / 3.0
    P_load_8_base = 50.0 / 3.0

    # Additional load (in MW)
    P_load_8_add = 100.0

    # Resistors (in Ohms)
    R = {
        'T7_8': 4.50,
        'T8_9': 6.30,
        'T5_7': 16.93,
        'T6_9': 20.63,
        'T4_5': 5.29,
        'T4_6': 8.99,
    }

    # Capacitor reactive power supply (in MVAr)
    Q_cap = 50.0

    # 2. Calculate Total Load
    P_load_total = P_load_5 + P_load_6 + P_load_8_base + P_load_8_add
    
    # Calculate reactive power from active power and power factor
    phi = math.acos(pf)
    tan_phi = math.tan(phi)
    Q_load_total = P_load_total * tan_phi

    # 3. Account for Reactive Power Compensation
    Q_net_demand = Q_load_total - Q_cap

    # 4. Determine Generator Dispatch (assuming equal sharing)
    # Total power to be generated
    S_gen_total = complex(P_load_total, Q_net_demand)

    # Power per generator (3 generators)
    S_gen_per = S_gen_total / 3.0
    
    # 5. Calculate Power Flow in Each Line
    
    # Complex power for each load
    S_load_5 = complex(P_load_5, P_load_5 * tan_phi)
    S_load_6 = complex(P_load_6, P_load_6 * tan_phi)
    P_load_8 = P_load_8_base + P_load_8_add
    S_load_8 = complex(P_load_8, P_load_8 * tan_phi)
    
    # Use symmetry: S_78 = S_98, S_75 = S_96, S_45 = S_46
    # Power balance at Bus 8: S_78 + S_98 = S_load_8
    S_78 = S_load_8 / 2.0
    S_98 = S_78
    
    # Power balance at Bus 7: S_gen_2 = S_75 + S_78
    # Generator 2 is connected to Bus 7
    S_75 = S_gen_per - S_78
    
    # Power balance at Bus 9: S_gen_3 = S_96 + S_98
    S_96 = S_gen_per - S_98 # This will be equal to S_75 by symmetry

    # Power balance at Bus 5: S_45 + S_75 = S_load_5
    S_45 = S_load_5 - S_75

    # Power balance at Bus 6: S_46 + S_96 = S_load_6
    S_46 = S_load_6 - S_96 # This will be equal to S_45 by symmetry

    # Power flow magnitudes squared (in MVA^2)
    S_mag_sq = {
        'T7_8': abs(S_78)**2,
        'T8_9': abs(S_98)**2,
        'T5_7': abs(S_75)**2,
        'T6_9': abs(S_96)**2,
        'T4_5': abs(S_45)**2,
        'T4_6': abs(S_46)**2,
    }

    # 6. Calculate Resistive Losses
    losses = {}
    V_sq = V_line_kV**2
    total_loss = 0
    
    print("Calculating individual line losses (in MW):")
    print("Formula: P_loss = |S|^2 * R / V^2\n")

    for line, S2 in S_mag_sq.items():
        R_val = R[line]
        loss = (S2 * R_val) / V_sq
        losses[line] = loss
        total_loss += loss
        print(f"Loss in {line}: (|{S2**0.5:.3f} MVA|^2 * {R_val:.2f} Ω) / ({V_line_kV:.1f} kV)^2 = {loss:.3f} MW")
    
    # 7. Sum Total Losses and Print Result
    print("\n-------------------------------------------")
    print(f"Total Power Loss = {' + '.join([f'{l:.3f}' for l in losses.values()])} MW")
    print(f"Total Power Loss = {total_loss:.3f} MW")
    
    # For automated checking
    return total_loss

if __name__ == '__main__':
    final_answer = calculate_power_losses()
    print(f"<<<{final_answer:.3f}>>>")