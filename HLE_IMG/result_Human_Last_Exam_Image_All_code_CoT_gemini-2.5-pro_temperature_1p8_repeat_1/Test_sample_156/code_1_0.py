import math

def calculate_reactive_power_deficit():
    """
    Calculates the reactive power deficit for the given power system.
    """

    # Step 1: Define initial values from the problem description
    
    # Initial Load Information
    P_load1 = 125.0  # MW
    Q_load1_initial = 50.0  # MVAr

    P_load2 = 90.0  # MW
    Q_load2_initial = 30.0  # MVAr

    P_load3 = 100.0  # MW
    Q_load3_initial = 35.0  # MVAr

    # Line Information
    line_reactance_per_km = 0.12  # Ohm/km
    line_length = 80.0  # km
    line_voltage = 230.0  # kV

    # Reactive Power Support
    s3_support_factor = 0.50  # 50%
    battery_support = 10.0  # MVAr

    # System condition factors
    reactive_demand_increase = 0.10  # 10%

    # Step 2: Calculate new reactive power demand from loads with the 10% increase
    q_increase_factor = 1 + reactive_demand_increase
    Q_load1_new = Q_load1_initial * q_increase_factor
    Q_load2_new = Q_load2_initial * q_increase_factor
    Q_load3_new = Q_load3_initial * q_increase_factor
    total_Q_demand_from_loads = Q_load1_new + Q_load2_new + Q_load3_new

    # Step 3: Calculate reactive power losses in the transmission line
    # Assumption: The power flowing on the Bus 7-9 line is to serve Load 3.
    total_line_reactance = line_reactance_per_km * line_length
    S_load3_new_squared = P_load3**2 + Q_load3_new**2
    Q_line_loss = total_line_reactance * S_load3_new_squared / (line_voltage**2)

    # Step 4: Calculate total system reactive demand
    total_Q_demand_system = total_Q_demand_from_loads + Q_line_loss

    # Step 5: Calculate total available reactive power supply
    # Assumption: S3 support is based on the initial load Q.
    Q_s3_support = s3_support_factor * Q_load1_initial
    total_Q_supply = Q_s3_support + battery_support

    # Step 6: Calculate the reactive power deficit
    deficit = total_Q_demand_system - total_Q_supply

    # Step 7: Print the detailed calculation
    print("--- Reactive Power Deficit Calculation ---\n")
    print("1. Total Reactive Power Demand Calculation:")
    print(f"  - New demand for Load 1 (10% increase): {Q_load1_new:.2f} MVAr")
    print(f"  - New demand for Load 2 (10% increase): {Q_load2_new:.2f} MVAr")
    print(f"  - New demand for Load 3 (10% increase): {Q_load3_new:.2f} MVAr")
    print(f"  - Line Reactive Power Loss (Q_loss): {Q_line_loss:.2f} MVAr")
    print(f"  - Total Reactive Demand = {Q_load1_new:.2f} + {Q_load2_new:.2f} + {Q_load3_new:.2f} + {Q_line_loss:.2f} = {total_Q_demand_system:.2f} MVAr\n")

    print("2. Total Reactive Power Supply Calculation:")
    print(f"  - Support from Generator S3 (50% of initial Q_Load1): {Q_s3_support:.2f} MVAr")
    print(f"  - Support from Battery: {battery_support:.2f} MVAr")
    print(f"  - Total Reactive Supply = {Q_s3_support:.2f} + {battery_support:.2f} = {total_Q_supply:.2f} MVAr\n")
    
    print("3. Final Deficit Calculation:")
    print(f"  - Deficit = Total Demand - Total Supply")
    print(f"  - Deficit = {total_Q_demand_system:.2f} MVAr - {total_Q_supply:.2f} MVAr")
    print(f"  - The total reactive power deficit is {deficit:.2f} MVAr.")

calculate_reactive_power_deficit()