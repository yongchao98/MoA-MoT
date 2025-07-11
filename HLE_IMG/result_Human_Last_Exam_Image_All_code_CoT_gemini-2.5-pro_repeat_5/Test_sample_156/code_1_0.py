import math

def calculate_power_deficit():
    """
    Calculates the reactive power deficit based on the given power system parameters.
    """
    # Step 1: Define initial values from the problem statement
    # Load Information (initial)
    P_load1 = 125  # MW
    Q_load1 = 50   # MVAr

    P_load2 = 90   # MW
    Q_load2 = 30   # MVAr

    P_load3 = 100  # MW
    Q_load3 = 35   # MVAr

    # Line parameters
    Z_per_km_imag = 0.12  # Ohm/km (Reactance)
    line_length = 80      # km
    line_voltage = 230    # kV

    # Reactive Power Support
    S3_support_factor = 0.50 # 50% of Load 1's reactive power
    BESS_support = 10        # MVAr

    # System changes
    q_demand_increase_factor = 1.10 # 10% increase

    # --- Calculations ---

    # Step 2: Calculate the new total reactive power demand from loads
    Q_load1_new = Q_load1 * q_demand_increase_factor
    Q_load2_new = Q_load2 * q_demand_increase_factor
    Q_load3_new = Q_load3 * q_demand_increase_factor
    Q_demand_total_loads = Q_load1_new + Q_load2_new + Q_load3_new

    # Step 3: Calculate the reactive power loss in the transmission line
    P_total = P_load1 + P_load2 + P_load3
    S_total_MVA = math.sqrt(P_total**2 + Q_demand_total_loads**2)
    X_line_total = Z_per_km_imag * line_length
    Q_loss_MVar = (S_total_MVA**2 * X_line_total) / (line_voltage**2)

    # Step 4: Calculate total reactive power need
    Q_total_need = Q_demand_total_loads + Q_loss_MVar

    # Step 5: Calculate total available reactive power from specified sources
    Q_support_S3 = Q_load1 * S3_support_factor
    Q_total_supply_specified = Q_support_S3 + BESS_support

    # Step 6: Calculate the reactive power deficit
    deficit_MVar = Q_total_need - Q_total_supply_specified

    # --- Print Results ---
    print("Calculation of Reactive Power Deficit")
    print("="*50)
    
    print("\n1. Total Updated Reactive Power Demand (Loads):")
    print(f"   Demand = (Q_Load1 * 1.1) + (Q_Load2 * 1.1) + (Q_Load3 * 1.1)")
    print(f"   Demand = ({Q_load1} * {q_demand_increase_factor}) + ({Q_load2} * {q_demand_increase_factor}) + ({Q_load3} * {q_demand_increase_factor}) = {Q_demand_total_loads:.2f} MVAr")

    print("\n2. Reactive Power Loss (Line):")
    print(f"   Total Real Power P = {P_load1} + {P_load2} + {P_load3} = {P_total} MW")
    print(f"   Total Apparent Power S = sqrt({P_total}^2 + {Q_demand_total_loads:.2f}^2) = {S_total_MVA:.2f} MVA")
    print(f"   Line Reactance X = {Z_per_km_imag} * {line_length} = {X_line_total:.2f} Ohm")
    print(f"   Loss = (S^2 * X) / V^2 = ({S_total_MVA:.2f}^2 * {X_line_total:.2f}) / {line_voltage}^2 = {Q_loss_MVar:.2f} MVAr")

    print("\n3. Total Specified Reactive Power Supply:")
    print(f"   Supply = (S3 Support) + (BESS Support)")
    print(f"   Supply = ({Q_load1} * {S3_support_factor}) + {BESS_support} = {Q_total_supply_specified:.2f} MVAr")

    print("\n" + "="*50)
    print("Final Deficit Calculation")
    print("="*50)
    print("Deficit = (Total Load Demand + Line Loss) - Specified Supply")
    print(f"Deficit = ({Q_demand_total_loads:.2f} MVAr + {Q_loss_MVar:.2f} MVAr) - {Q_total_supply_specified:.2f} MVAr")
    print(f"Deficit = {Q_total_need:.2f} MVAr - {Q_total_supply_specified:.2f} MVAr = {deficit_MVar:.2f} MVAr")

calculate_power_deficit()