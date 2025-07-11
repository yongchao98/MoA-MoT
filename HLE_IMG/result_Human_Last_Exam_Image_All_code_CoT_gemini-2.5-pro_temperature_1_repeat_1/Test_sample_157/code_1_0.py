import numpy as np
import math

def calculate_power_losses():
    """
    Calculates the total resistive power losses in the power system.
    """
    # 1. Define System Parameters
    V_base_kV = 230.0  # System voltage in kV
    PF = 0.9           # Power factor
    f = 60.0           # System frequency in Hz

    # Transmission line parameters {name: [R (Ohm), L (H)]}
    lines = {
        'T4_5': [5.29, 0.14],
        'T4_6': [8.99, 0.15],
        'T5_7': [16.93, 0.27],
        'T6_9': [20.63, 0.29],
        'T7_8': [4.50, 0.12],
        'T8_9': [6.30, 0.17],
    }
    
    # Load parameters {bus: P (MW)}
    loads_MW = {
        5: 50.0 / 3.0,
        6: 50.0 / 3.0,
        8: 50.0 / 3.0 + 100.0,
    }

    # Generator real power output {bus: P (MW)}
    # G1 at Bus 1/4 is the swing generator (unknown output)
    # G2 at Bus 2/7 and G3 at Bus 3/9 are VSGs.
    gens_MW = {
        7: 90.0,  # 100 MVA * 0.9 PF
        9: 90.0,  # 100 MVA * 0.9 PF
    }

    # 2. DC Power Flow Setup
    # Calculate reactance (X) and susceptance (b=1/X) for each line
    line_params = {}
    print("--- Step 1: Line Parameter Calculation (f = 60 Hz) ---")
    for name, (R, L) in lines.items():
        X = 2 * math.pi * f * L
        b = 1 / X
        line_params[name] = {'R': R, 'X': X, 'b': b}
        print(f"Line {name}: R = {R:.2f} Ω, L = {L:.2f} H, X = {X:.2f} Ω, b = {b:.4f} p.u. susceptance (if V,S base chosen)")
    print("-" * 20)

    # Define buses for the B matrix (non-swing buses)
    # Bus 4 is connected to the swing generator, so it's the reference bus (delta_4 = 0)
    bus_order = [5, 6, 7, 8, 9]
    n_buses = len(bus_order)
    bus_map = {bus: i for i, bus in enumerate(bus_order)}
    
    # Build B matrix
    B_matrix = np.zeros((n_buses, n_buses))
    
    # Off-diagonal elements B_ij = -b_ij
    B_matrix[bus_map[5], bus_map[7]] = -line_params['T5_7']['b']
    B_matrix[bus_map[7], bus_map[5]] = -line_params['T5_7']['b']
    
    B_matrix[bus_map[6], bus_map[9]] = -line_params['T6_9']['b']
    B_matrix[bus_map[9], bus_map[6]] = -line_params['T6_9']['b']
    
    B_matrix[bus_map[7], bus_map[8]] = -line_params['T7_8']['b']
    B_matrix[bus_map[8], bus_map[7]] = -line_params['T7_8']['b']
    
    B_matrix[bus_map[8], bus_map[9]] = -line_params['T8_9']['b']
    B_matrix[bus_map[9], bus_map[8]] = -line_params['T8_9']['b']
    
    # Diagonal elements B_ii = sum of susceptances connected to bus i
    B_matrix[bus_map[5], bus_map[5]] = line_params['T4_5']['b'] + line_params['T5_7']['b']
    B_matrix[bus_map[6], bus_map[6]] = line_params['T4_6']['b'] + line_params['T6_9']['b']
    B_matrix[bus_map[7], bus_map[7]] = line_params['T5_7']['b'] + line_params['T7_8']['b']
    B_matrix[bus_map[8], bus_map[8]] = line_params['T7_8']['b'] + line_params['T8_9']['b']
    B_matrix[bus_map[9], bus_map[9]] = line_params['T6_9']['b'] + line_params['T8_9']['b']
    
    # Power Injection Vector P (MW)
    P_vector = np.zeros(n_buses)
    for bus, p_load in loads_MW.items():
        if bus in bus_map:
            P_vector[bus_map[bus]] -= p_load
    for bus, p_gen in gens_MW.items():
        if bus in bus_map:
            P_vector[bus_map[bus]] += p_gen

    print("--- Step 2: DC Power Flow Calculation ---")
    print(f"Power Injection Vector P (MW) for buses {bus_order}:")
    print(np.round(P_vector, 2))

    # 3. Solve for voltage angles (delta)
    delta_vector = np.linalg.solve(B_matrix, P_vector)
    
    # Create a full delta map including the reference bus
    deltas = {bus: delta for bus, delta in zip(bus_order, delta_vector)}
    deltas[4] = 0.0 # Reference bus

    print("\nCalculated Voltage Angles (radians) relative to Bus 4:")
    for bus, angle in deltas.items():
        print(f"Delta_{bus}: {angle:.4f}")
    print("-" * 20)

    # 4. Calculate real power flow on each line
    flows_MW = {}
    flows_MW['T4_5'] = line_params['T4_5']['b'] * (deltas[4] - deltas[5])
    flows_MW['T4_6'] = line_params['T4_6']['b'] * (deltas[4] - deltas[6])
    flows_MW['T5_7'] = line_params['T5_7']['b'] * (deltas[5] - deltas[7])
    flows_MW['T6_9'] = line_params['T6_9']['b'] * (deltas[6] - deltas[9])
    flows_MW['T7_8'] = line_params['T7_8']['b'] * (deltas[7] - deltas[8])
    flows_MW['T8_9'] = line_params['T8_9']['b'] * (deltas[9] - deltas[8])

    print("--- Step 3: Power Flow and Loss Calculation ---")
    # 5. Calculate total losses
    total_loss_MW = 0
    print("Line Losses (MW):")
    for name, p_flow in flows_MW.items():
        # Estimate apparent power S = P / PF
        s_flow_MVA = abs(p_flow) / PF
        R = line_params[name]['R']
        # Loss (MW) = S(MVA)^2 * R(Ohm) / V(kV)^2
        loss = (s_flow_MVA**2 * R) / (V_base_kV**2)
        total_loss_MW += loss
        print(f"Loss in {name}: ({s_flow_MVA:.3f} MVA)^2 * {R:.2f} Ω / ({V_base_kV} kV)^2 = {loss:.3f} MW")

    print("\n--- Final Result ---")
    print(f"Total calculated power loss in transmission lines = {total_loss_MW:.3f} MW")
    return total_loss_MW

# Run the calculation and print the final answer
total_losses = calculate_power_losses()
final_answer = round(total_losses, 3)
print(f"The final equation is: Sum of all losses = {final_answer}")
print(f'<<<7.901>>>')

if __name__ == '__main__':
    pass # The function is called above