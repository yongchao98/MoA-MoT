import math

def calculate_reactive_power_deficit():
    """
    Calculates the total reactive power deficit for the given power system.
    The function follows these steps:
    1.  Calculates the total reactive power demand from loads, including a 10% increase.
    2.  Calculates the reactive power losses in the transmission line from Bus 7 to Bus 9, assuming it's tapped at the midpoint (Bus 8).
    3.  Sums the load demand and line losses to get the total system demand.
    4.  Calculates the total available reactive power from specified sources (Generator S3 and BESS).
    5.  Calculates the deficit by subtracting the total supply from the total demand.
    """
    # Step 1: Define initial values from the problem description and diagram
    # Load information (initial values)
    P_load1 = 125.0  # MW
    Q_load1_initial = 50.0  # MVAr

    P_load2 = 90.0   # MW
    Q_load2_initial = 30.0  # MVAr

    P_load3 = 100.0  # MW
    Q_load3_initial = 35.0  # MVAr

    # Line information for the line between Bus 7 and Bus 9
    line_impedance_imag_per_km = 0.12  # Ohm/km
    total_line_length = 80.0  # km
    line_voltage_kv = 230.0  # kV

    # System factors
    voltage_drop_q_increase = 0.10  # 10% increase in reactive power demand

    # Reactive power support
    s3_support_factor = 0.50  # Generator S3 provides 50% of Load 1's reactive power
    bess_support = 10.0  # MVAr from battery

    # Assumption: Bus 8 is at the midpoint of the line between Bus 7 and Bus 9.
    length_segment_1 = total_line_length / 2.0 # km
    length_segment_2 = total_line_length / 2.0 # km

    # Step 2: Calculate the new reactive power demand for each load
    q_increase_multiplier = 1 + voltage_drop_q_increase
    Q_load1_final = Q_load1_initial * q_increase_multiplier
    Q_load2_final = Q_load2_initial * q_increase_multiplier
    Q_load3_final = Q_load3_initial * q_increase_multiplier

    total_Q_loads_final = Q_load1_final + Q_load2_final + Q_load3_final

    # Step 3: Calculate reactive power losses in the transmission line segments
    # Reactance of each segment
    X_segment_1 = line_impedance_imag_per_km * length_segment_1
    X_segment_2 = line_impedance_imag_per_km * length_segment_2

    # Power flow and loss in segment 1 (Bus 7 -> Bus 8), carrying power for Loads 2 & 3
    P_flow_1 = P_load2 + P_load3
    Q_flow_1 = Q_load2_final + Q_load3_final
    S_flow_1_mva = math.sqrt(P_flow_1**2 + Q_flow_1**2)
    I_1_amps = (S_flow_1_mva * 1e6) / (line_voltage_kv * 1e3)
    Q_loss_1_mvar = (I_1_amps**2 * X_segment_1) / 1e6

    # Power flow and loss in segment 2 (Bus 8 -> Bus 9), carrying power for Load 2
    P_flow_2 = P_load2
    Q_flow_2 = Q_load2_final
    S_flow_2_mva = math.sqrt(P_flow_2**2 + Q_flow_2**2)
    I_2_amps = (S_flow_2_mva * 1e6) / (line_voltage_kv * 1e3)
    Q_loss_2_mvar = (I_2_amps**2 * X_segment_2) / 1e6

    # Total reactive power loss
    total_Q_loss = Q_loss_1_mvar + Q_loss_2_mvar

    # Step 4: Calculate total reactive power demand
    total_Q_demand = total_Q_loads_final + total_Q_loss

    # Step 5: Calculate total available reactive power supply
    Q_from_S3 = s3_support_factor * Q_load1_initial
    total_Q_supply = Q_from_S3 + bess_support

    # Step 6: Calculate the reactive power deficit
    deficit = total_Q_demand - total_Q_supply

    # Step 7: Print the final calculation
    print("Calculation of Reactive Power Deficit")
    print("---------------------------------------")
    print(f"Total Reactive Power Demand = (Increased Load Demand) + (Line Losses)")
    print(f"Total Reactive Power Demand = ({total_Q_loads_final:.2f} MVAr) + ({total_Q_loss:.2f} MVAr) = {total_Q_demand:.2f} MVAr\n")
    
    print(f"Total Available Supply = (Supply from S3) + (Supply from BESS)")
    print(f"Total Available Supply = ({Q_from_S3:.2f} MVAr) + ({bess_support:.2f} MVAr) = {total_Q_supply:.2f} MVAr\n")

    print("Final Equation:")
    print(f"Reactive Power Deficit = Total Demand - Total Supply")
    print(f"Reactive Power Deficit = {total_Q_demand:.2f} - {total_Q_supply:.2f} = {deficit:.2f} MVAr")

if __name__ == "__main__":
    calculate_reactive_power_deficit()
    # The final numerical answer for the deficit is 96.07
    print("\n<<<96.07>>>")