import math

def calculate_reactive_power_deficit():
    """
    Calculates the total reactive power deficit based on the given power system parameters.
    """
    # --- Step 1: Define Initial Parameters ---
    # Load Information (in MW and MVAr)
    P_L1 = 125.0
    Q_L1_initial = 50.0
    P_L2 = 90.0
    Q_L2_initial = 30.0
    P_L3 = 100.0
    Q_L3_initial = 35.0

    # Line Information
    Z_line_per_km_imag = 0.12  # Ohm/km (reactance)
    line_length = 80.0  # km
    line_voltage = 230.0  # kV

    # System Factors
    voltage_drop_q_increase = 0.10  # 10%

    # Reactive Power Support
    S3_support_factor = 0.50  # 50%
    BESS_support = 10.0  # MVAr

    print("Step-by-step calculation for reactive power deficit:")
    print("-" * 50)

    # --- Step 2: Calculate Final Reactive Power Demand from Loads ---
    print("1. Calculate final reactive power demand from loads (with 10% increase):")
    Q_L1_final = Q_L1_initial * (1 + voltage_drop_q_increase)
    Q_L2_final = Q_L2_initial * (1 + voltage_drop_q_increase)
    Q_L3_final = Q_L3_initial * (1 + voltage_drop_q_increase)
    print(f"  - Final Q for Load 1: {Q_L1_initial:.1f} MVAr * 1.10 = {Q_L1_final:.1f} MVAr")
    print(f"  - Final Q for Load 2: {Q_L2_initial:.1f} MVAr * 1.10 = {Q_L2_final:.1f} MVAr")
    print(f"  - Final Q for Load 3: {Q_L3_initial:.1f} MVAr * 1.10 = {Q_L3_final:.1f} MVAr")

    total_load_q = Q_L1_final + Q_L2_final + Q_L3_final
    print(f"  - Total Reactive Power Demand from Loads: {Q_L1_final:.1f} + {Q_L2_final:.1f} + {Q_L3_final:.1f} = {total_load_q:.1f} MVAr\n")

    # --- Step 3: Calculate Reactive Power Loss in the Transmission Line ---
    print("2. Calculate reactive power loss in the transmission line:")
    # Power assumed to flow to loads 2 and 3
    P_23 = P_L2 + P_L3
    Q_23 = Q_L2_final + Q_L3_final
    
    # Apparent power S squared for these loads
    S_23_sq = P_23**2 + Q_23**2
    S_23 = math.sqrt(S_23_sq)
    print(f"  - Apparent Power to Loads 2&3 (S) = sqrt(({P_L2:.1f} + {P_L3:.1f})^2 + ({Q_L2_final:.1f} + {Q_L3_final:.1f})^2) = {S_23:.2f} MVA")

    # Total line reactance
    X_line = Z_line_per_km_imag * line_length
    print(f"  - Total Line Reactance (X) = {Z_line_per_km_imag:.2f} Ohm/km * {line_length:.1f} km = {X_line:.2f} Ohms")

    # Reactive power loss Q_loss = |S|^2 * X / |V|^2
    Q_loss = S_23_sq * X_line / (line_voltage**2)
    print(f"  - Reactive Power Loss (Q_loss) = ({S_23:.2f} MVA)^2 * {X_line:.2f} Ohms / ({line_voltage:.1f} kV)^2 = {Q_loss:.2f} MVAr\n")

    # --- Step 4: Calculate Total Reactive Power Demand ---
    print("3. Calculate total reactive power demand (Loads + Line Loss):")
    total_q_demand = total_load_q + Q_loss
    print(f"  - Total Demand = {total_load_q:.1f} MVAr (Loads) + {Q_loss:.2f} MVAr (Loss) = {total_q_demand:.2f} MVAr\n")

    # --- Step 5: Calculate Total Available Reactive Power Supply ---
    print("4. Calculate total available reactive power supply:")
    q_support_s3 = S3_support_factor * Q_L1_initial
    print(f"  - S3 Support = {S3_support_factor:.2f} * {Q_L1_initial:.1f} MVAr = {q_support_s3:.1f} MVAr")
    print(f"  - BESS Support = {BESS_support:.1f} MVAr")

    total_q_supply = q_support_s3 + BESS_support
    print(f"  - Total Supply = {q_support_s3:.1f} MVAr + {BESS_support:.1f} MVAr = {total_q_supply:.1f} MVAr\n")

    # --- Step 6: Calculate the Reactive Power Deficit ---
    print("5. Calculate the final reactive power deficit:")
    deficit = total_q_demand - total_q_supply
    print(f"  - Deficit = Total Demand - Total Supply")
    
    # Final Equation Output
    print("\nFinal Equation:")
    print(f"Deficit = (({Q_L1_final:.1f}) + ({Q_L2_final:.1f}) + ({Q_L3_final:.1f}) + ({Q_loss:.2f})) - (({q_support_s3:.1f}) + ({BESS_support:.1f}))")
    print(f"Deficit = ({total_q_demand:.2f}) - ({total_q_supply:.1f})")
    print(f"Deficit = {deficit:.2f} MVAr")
    
    return deficit

if __name__ == '__main__':
    final_deficit = calculate_reactive_power_deficit()
    # The final answer is printed within the function.
    # For the required format, we print it again here.
    # print(f"<<<{final_deficit:.2f}>>>")
    
calculate_reactive_power_deficit()
<<<98.98>>>