import cmath
import math

def calculate_power_losses():
    """
    Calculates the total resistive power losses in the power system.
    """
    # Step 1: Define System Parameters
    # System frequency (assuming 60 Hz as it's standard and not given)
    f = 60
    omega = 2 * math.pi * f

    # Line Voltages
    V_ll_hv = 230e3  # Line-to-line voltage on the high voltage side (V)

    # Power Factor
    pf = 0.9
    angle = math.acos(pf)
    tan_phi = math.tan(angle)

    # Transmission Line Impedances (R in Ohm, L in Henry)
    lines = {
        'T7_8': {'R': 4.50, 'L': 0.12},
        'T8_9': {'R': 6.30, 'L': 0.17},
        'T5_7': {'R': 16.93, 'L': 0.27},
        'T6_9': {'R': 20.63, 'L': 0.29},
        'T4_5': {'R': 5.29, 'L': 0.14},
        'T4_6': {'R': 8.99, 'L': 0.15},
    }
    # Calculate complex impedance Z = R + j*omega*L for each line
    Z = {name: data['R'] + 1j * omega * data['L'] for name, data in lines.items()}

    # Step 2: Calculate Complex Power for Loads
    # Active Power Loads (MW)
    P_L5 = 50/3 * 1e6
    P_L6 = 50/3 * 1e6
    P_L8_base = 50/3 * 1e6
    P_L8_add = 100 * 1e6
    P_L8 = P_L8_base + P_L8_add

    # Complex Power Loads (S = P + jQ)
    S_L5 = P_L5 + 1j * P_L5 * tan_phi
    S_L6 = P_L6 + 1j * P_L6 * tan_phi
    S_L8 = P_L8 + 1j * P_L8 * tan_phi

    # Capacitor Bank supplies 50 MVAr
    S_cap = -1j * 50e6

    # Step 3: Power Distribution Assumption
    # Assume generators share the net load equally.
    S_load_total = S_L5 + S_L6 + S_L8
    S_net_demand = S_load_total + S_cap # S_cap is negative (supply)
    S_gen_per_unit = S_net_demand / 3

    # Power injected into the grid from each generator complex
    S_inj_7 = S_gen_per_unit  # From VSG 2 at Bus 7
    S_inj_9 = S_gen_per_unit  # From VSG 3 at Bus 9
    # From Gen 1 + Capacitor at Bus 4
    S_inj_4_eff = S_gen_per_unit - S_cap # Gen 1 compensates for cap being at its bus

    # Step 4: Solve Power Network for Line Flows
    # The system of equations was solved as follows:
    # Power to Bus 8 is split between lines T7_8 and T8_9
    # based on their relative impedances.
    Z_sum_78_89 = Z['T7_8'] + Z['T8_9']
    S_78 = S_L8 * Z['T8_9'] / Z_sum_78_89
    S_98 = S_L8 * Z['T7_8'] / Z_sum_78_89 # Flow from 9 to 8

    # Power flows out of buses 7 and 9
    S_75 = S_inj_7 - S_78
    S_96 = S_inj_9 - S_98

    # Power flows out of bus 4
    S_45 = S_L5 - S_75
    S_46 = S_L6 - S_96

    # Store all line flows
    S_flows = {
        'T7_8': S_78,
        'T8_9': -S_98, # Flow from 8 to 9
        'T5_7': -S_75, # Flow from 5 to 7
        'T6_9': -S_96, # Flow from 6 to 9
        'T4_5': S_45,
        'T4_6': S_46
    }

    # Step 5: Calculate Resistive Losses for Each Line
    losses = {}
    V_sq = V_ll_hv**2
    for name, s_flow in S_flows.items():
        R_line = lines[name]['R']
        loss_watts = (abs(s_flow)**2 / V_sq) * R_line
        losses[name] = loss_watts / 1e6 # Convert to MW

    # Step 6: Summation and Final Output
    total_loss = sum(losses.values())

    print("Individual resistive power losses:")
    loss_strings = []
    for name, loss in losses.items():
        print(f"Loss in {name}: {loss:.3f} MW")
        loss_strings.append(f"{loss:.3f}")
        
    print("\nFinal Calculation:")
    print(f"Total Loss = {' + '.join(loss_strings)} = {total_loss:.3f} MW")

calculate_power_losses()
<<<1.622>>>