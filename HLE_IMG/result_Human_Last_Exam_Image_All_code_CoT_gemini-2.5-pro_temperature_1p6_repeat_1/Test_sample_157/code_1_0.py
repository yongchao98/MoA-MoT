import math

def calculate_power_loss():
    """
    Calculates the total power loss in the system based on generation and load.
    """

    # --- System Generation Parameters ---
    # Three generators, each rated at 100 MVA
    s_gen_1 = 100.0  # MVA
    s_gen_2 = 100.0  # MVA
    s_gen_3 = 100.0  # MVA

    # System-wide power factor
    power_factor = 0.9

    # --- System Load Parameters ---
    # Three base loads, each 50/3 MW
    p_base_load = 50.0 / 3.0  # MW
    num_base_loads = 3

    # Additional load is active, 100 MW
    p_additional_load = 100.0  # MW

    # --- Step 1: Calculate Total Generated Active Power ---
    # Total apparent power generation at 100% capacity
    s_gen_total = s_gen_1 + s_gen_2 + s_gen_3
    # Total active power generation
    p_gen_total = s_gen_total * power_factor

    print("--- Power Generation Calculation ---")
    print(f"Total Apparent Power Generated (S_gen) = {s_gen_1:.1f} MVA + {s_gen_2:.1f} MVA + {s_gen_3:.1f} MVA = {s_gen_total:.3f} MVA")
    print(f"Total Active Power Generated (P_gen) = S_gen * PF = {s_gen_total:.3f} * {power_factor} = {p_gen_total:.3f} MW\n")

    # --- Step 2: Calculate Total Load Active Power ---
    p_base_total = p_base_load * num_base_loads
    p_load_total = p_base_total + p_additional_load

    print("--- Power Load Calculation ---")
    print(f"Total Base Load = 3 * ({p_base_load:.3f} MW) = {p_base_total:.3f} MW")
    print(f"Additional Active Load = {p_additional_load:.3f} MW")
    print(f"Total Active Power Load (P_load) = {p_base_total:.3f} MW + {p_additional_load:.3f} MW = {p_load_total:.3f} MW\n")

    # --- Step 3: Calculate Total Power Loss ---
    # Power loss is the difference between generation and load
    p_loss = p_gen_total - p_load_total

    print("--- Power Loss Calculation ---")
    print("The total power loss is the difference between total generated power and total consumed power.")
    print(f"Power Loss = P_gen_total - P_load_total")
    print(f"Power Loss = {p_gen_total:.3f} MW - {p_load_total:.3f} MW = {p_loss:.3f} MW")
    
    # Final answer
    print(f"\nThe final calculated total power loss is {p_loss:.3f} MW.")

# Execute the calculation
calculate_power_loss()