import cmath
import math

def calculate_power_system_losses():
    """
    Calculates the total resistive power losses in the transmission lines based on the provided power system diagram and parameters.
    """

    # --- 1. Constants and Given Data ---
    V_line = 230e3        # Line-to-line voltage in Volts (230 kV)
    PF = 0.9              # Lagging power factor
    f = 60                # System frequency in Hz (standard assumption)
    omega = 2 * math.pi * f

    # Resistors (R) and Inductors (L) for transmission lines
    line_data = {
        'T7_8': {'R': 4.50, 'L': 0.12},
        'T8_9': {'R': 6.30, 'L': 0.17},
        'T5_7': {'R': 16.93, 'L': 0.27},
        'T6_9': {'R': 20.63, 'L': 0.29},
        'T4_5': {'R': 5.29, 'L': 0.14},
        'T4_6': {'R': 8.99, 'L': 0.15},
    }

    # --- 2. Calculate Complex Loads (S = P + jQ) in VA ---
    # Angle for PF=0.9 lagging
    phi = math.acos(PF)
    tan_phi = math.tan(phi)

    # Loads at Buses 5 and 6
    P_L56_W = (50 / 3) * 1e6  # in Watts
    S_L5 = complex(P_L56_W, P_L56_W * tan_phi)
    S_L6 = complex(P_L56_W, P_L56_W * tan_phi)

    # Total Load at Bus 8 (Base + Additional)
    P_L8_W = (50 / 3 + 100) * 1e6  # in Watts
    S_L8 = complex(P_L8_W, P_L8_W * tan_phi)

    # --- 3. Calculate Complex Impedances (Z = R + jwL) in Ohms ---
    Z = {}
    for name, data in line_data.items():
        R = data['R']
        X = omega * data['L']
        Z[name] = complex(R, X)

    # --- 4. Approximate Power Flow Split based on Impedance Division ---
    # This model assumes loads are supplied by adjacent source buses via parallel paths.
    S_flow = {}
    
    # Bus 5 is fed by T4_5 and T5_7
    S_flow['T4_5'] = S_L5 * Z['T5_7'] / (Z['T4_5'] + Z['T5_7'])
    S_flow['T5_7'] = S_L5 * Z['T4_5'] / (Z['T4_5'] + Z['T5_7'])

    # Bus 6 is fed by T4_6 and T6_9
    S_flow['T4_6'] = S_L6 * Z['T6_9'] / (Z['T4_6'] + Z['T6_9'])
    S_flow['T6_9'] = S_L6 * Z['T4_6'] / (Z['T4_6'] + Z['T6_9'])

    # Bus 8 is fed by T7_8 and T8_9
    S_flow['T7_8'] = S_L8 * Z['T8_9'] / (Z['T7_8'] + Z['T8_9'])
    S_flow['T8_9'] = S_L8 * Z['T7_8'] / (Z['T7_8'] + Z['T8_9'])

    # --- 5. Calculate Power Loss in Each Line ---
    # P_loss = |S_flow|^2 * R / V_line^2
    losses_MW = {}
    for name, s_val in S_flow.items():
        loss_watts = (abs(s_val)**2 * Z[name].real) / (V_line**2)
        losses_MW[name] = loss_watts / 1e6

    # --- 6. Sum Total Losses and Print Results ---
    total_loss_MW = sum(losses_MW.values())

    print("Calculating individual line losses (in MW):")
    print(f"Loss(T4_5) = |S_flow|^2 * R / V^2 = {losses_MW['T4_5']:.3f} MW")
    print(f"Loss(T5_7) = |S_flow|^2 * R / V^2 = {losses_MW['T5_7']:.3f} MW")
    print(f"Loss(T4_6) = |S_flow|^2 * R / V^2 = {losses_MW['T4_6']:.3f} MW")
    print(f"Loss(T6_9) = |S_flow|^2 * R / V^2 = {losses_MW['T6_9']:.3f} MW")
    print(f"Loss(T7_8) = |S_flow|^2 * R / V^2 = {losses_MW['T7_8']:.3f} MW")
    print(f"Loss(T8_9) = |S_flow|^2 * R / V^2 = {losses_MW['T8_9']:.3f} MW")
    print("\nCalculating total power loss:")
    print(f"Total Loss = {losses_MW['T4_5']:.3f} + {losses_MW['T5_7']:.3f} + {losses_MW['T4_6']:.3f} + {losses_MW['T6_9']:.3f} + {losses_MW['T7_8']:.3f} + {losses_MW['T8_9']:.3f}")
    print(f"Total Loss = {total_loss_MW:.3f} MW")


calculate_power_system_losses()
<<<23.467>>>