import math

def calculate_power_loss():
    """
    Calculates the total resistive power losses in the power system based on the provided parameters.
    """

    # --- Given Parameters ---

    # Generation
    S_gen_1 = 100  # MVA, Synchronous Generator 1
    S_gen_2 = 100  # MVA, VSG Inverter 2
    S_gen_3 = 100  # MVA, VSG Inverter 3
    PF = 0.9      # Power factor for all components, lagging

    # Loads
    P_load_base = 50 / 3  # MW, for each of the 3 base loads
    num_base_loads = 3
    P_load_additional = 100  # MW, active additional load at Bus 8

    # --- Step 1: Calculate Total Generated Active Power ---
    
    print("Step 1: Calculate Total Generated Active Power (P_G)")
    
    # Total generated apparent power S_G is the sum of the capacities, as the system runs at 100%
    S_gen_total = S_gen_1 + S_gen_2 + S_gen_3
    print(f"The system operates at 100% capacity.")
    print(f"Total apparent power generation S_G = {S_gen_1} MVA (Gen 1) + {S_gen_2} MVA (VSG 2) + {S_gen_3} MVA (VSG 3) = {S_gen_total} MVA")
    
    # Total generated active power P_G is S_G multiplied by the power factor
    P_gen_total = S_gen_total * PF
    print(f"Given a system power factor of {PF} lagging.")
    print(f"Total active power generation P_G = S_G * PF = {S_gen_total} MVA * {PF} = {P_gen_total:.3f} MW\n")

    # --- Step 2: Calculate Total Consumed Active Power (Load) ---
    
    print("Step 2: Calculate Total Consumed Active Power (P_L)")
    
    # Total base load power
    P_load_base_total = num_base_loads * P_load_base
    
    # Total load is the sum of all base loads and the active additional load
    P_load_total = P_load_base_total + P_load_additional
    print(f"Total base load = {num_base_loads} * {P_load_base:.3f} MW = {P_load_base_total:.3f} MW")
    print(f"Additional active load at Bus 8 = {P_load_additional:.3f} MW")
    print(f"Total active power load P_L = {P_load_base_total:.3f} MW + {P_load_additional:.3f} MW = {P_load_total:.3f} MW\n")

    # --- Step 3: Calculate Total Power Loss ---
    
    print("Step 3: Calculate Total Resistive Power Loss (P_loss)")
    
    # By the principle of power balance, Total Power Loss = Total Generated Power - Total Load Power
    P_loss_total = P_gen_total - P_load_total
    print("The total power loss in the system is the difference between the total active power generated and the total active power consumed by loads.")
    print(f"P_loss = P_G - P_L")
    print(f"P_loss = {P_gen_total:.3f} MW - {P_load_total:.3f} MW = {P_loss_total:.3f} MW")
    
    # --- Final Answer ---
    
    final_answer = round(P_loss_total, 3)
    return final_answer

if __name__ == '__main__':
    loss = calculate_power_loss()
    print(f"<<<{loss:.3f}>>>")
