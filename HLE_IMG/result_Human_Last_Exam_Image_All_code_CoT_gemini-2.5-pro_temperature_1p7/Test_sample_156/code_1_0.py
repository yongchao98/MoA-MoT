import math

def calculate_power_deficit():
    """
    Calculates the reactive power deficit based on the provided power system parameters.
    """

    # --- Step 1: Define Given Data ---
    # Line Parameters
    Z_line_per_km_imag = 0.12  # Ohm/km
    line_length = 80  # km
    V_line_kV = 230  # kV

    # Initial Load Information
    P_L1, Q_L1_initial = 125, 50  # MW, MVAr
    P_L2, Q_L2_initial = 90, 30   # MW, MVAr
    P_L3, Q_L3_initial = 100, 35  # MW, MVAr

    # Support and System Effects
    S3_support_factor = 0.50
    BESS_support_mvar = 10
    q_demand_increase_factor = 0.10

    # --- Step 2: Calculate New Total Load Demand ---
    q_l1_new = Q_L1_initial * (1 + q_demand_increase_factor)
    q_l2_new = Q_L2_initial * (1 + q_demand_increase_factor)
    q_l3_new = Q_L3_initial * (1 + q_demand_increase_factor)
    total_q_demand_from_loads = q_l1_new + q_l2_new + q_l3_new

    # --- Step 3: Calculate Reactive Power Loss in the Line ---
    # a) Determine power that must flow from Bus 7 to Bus 9 area
    # This is the demand from loads L2 & L3 minus the local supply from S3 & BESS
    p_flow = P_L2 + P_L3
    q_supply_at_bus9 = (S3_support_factor * Q_L1_initial) + BESS_support_mvar
    q_demand_at_bus9_area = Q_L2_initial + Q_L3_initial
    q_flow = q_demand_at_bus9_area - q_supply_at_bus9
    
    # b) Calculate apparent power flow and current
    s_flow_mva = math.sqrt(p_flow**2 + q_flow**2)
    current_amps = (s_flow_mva * 1e6) / (math.sqrt(3) * V_line_kV * 1e3)

    # c) Calculate total line reactance
    total_line_reactance_ohm = Z_line_per_km_imag * line_length

    # d) Calculate 3-phase reactive power loss (Q_loss = 3 * I^2 * X)
    q_loss_var = 3 * (current_amps**2) * total_line_reactance_ohm
    q_loss_mvar = q_loss_var / 1e6

    # --- Step 4: Calculate Total Requirement and Supply ---
    total_q_requirement = total_q_demand_from_loads + q_loss_mvar
    total_q_supply = q_supply_at_bus9 # These are the only specified sources

    # --- Step 5: Calculate Deficit ---
    deficit = total_q_requirement - total_q_supply
    
    # --- Final Output ---
    print("Calculation of the Reactive Power Deficit:")
    print("1. Total Increased Reactive Load Demand:")
    print(f"   Load 1: {Q_L1_initial} MVAr * (1 + {q_demand_increase_factor}) = {q_l1_new:.1f} MVAr")
    print(f"   Load 2: {Q_L2_initial} MVAr * (1 + {q_demand_increase_factor}) = {q_l2_new:.1f} MVAr")
    print(f"   Load 3: {Q_L3_initial} MVAr * (1 + {q_demand_increase_factor}) = {q_l3_new:.1f} MVAr")
    
    print("\n2. Reactive Power Line Loss:")
    print(f"   Net Power Flow for Loss Calculation: {p_flow} MW, {q_flow} MVAr")
    print(f"   Calculated Line Loss (I²X): {q_loss_mvar:.2f} MVAr")
    
    print("\n3. Total Available Reactive Supply:")
    q_s3_supply = S3_support_factor * Q_L1_initial
    print(f"   From Generator S3: {S3_support_factor} * {Q_L1_initial} MVAr = {q_s3_supply:.1f} MVAr")
    print(f"   From BESS: {BESS_support_mvar:.1f} MVAr")
    
    print("\n--- Final Equation ---")
    # Deficit = (Total Increased Load Demand + Line Loss) - Total Available Supply
    print(f"Deficit = ({q_l1_new:.1f} + {q_l2_new:.1f} + {q_l3_new:.1f} + {q_loss_mvar:.2f}) - ({q_s3_supply:.1f} + {BESS_support_mvar:.1f})")
    print(f"Deficit = ({total_q_requirement:.2f}) - ({total_q_supply:.1f}) = {deficit:.2f} MVAr")

if __name__ == '__main__':
    calculate_power_deficit()
    # To directly output the final number for the answer block.
    # Recalculating here to isolate the final value.
    q_l1_new = 50 * 1.1
    q_l2_new = 30 * 1.1
    q_l3_new = 35 * 1.1
    total_q_demand_from_loads = q_l1_new + q_l2_new + q_l3_new
    p_flow = 90 + 100
    q_s3_supply = 0.5 * 50
    q_supply_at_bus9 = q_s3_supply + 10
    q_demand_at_bus9_area = 30 + 35
    q_flow = q_demand_at_bus9_area - q_supply_at_bus9
    s_flow_mva = math.sqrt(p_flow**2 + q_flow**2)
    current_amps = (s_flow_mva * 1e6) / (math.sqrt(3) * 230 * 1e3)
    total_line_reactance_ohm = 0.12 * 80
    q_loss_mvar = (3 * (current_amps**2) * total_line_reactance_ohm) / 1e6
    total_q_requirement = total_q_demand_from_loads + q_loss_mvar
    deficit = total_q_requirement - q_supply_at_bus9
    # The final print statement for the <<<answer>>> format is just the numerical result.
    final_answer = f"<<<{deficit:.2f}>>>"
