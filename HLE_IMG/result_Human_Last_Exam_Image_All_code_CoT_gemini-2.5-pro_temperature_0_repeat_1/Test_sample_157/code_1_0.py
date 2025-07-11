import math

def calculate_power_losses():
    """
    Calculates the total power losses in the system based on the provided data.
    """
    # --- Given Parameters ---
    # Generator parameters
    s_gen_per_unit = 100.0  # MVA per generator
    num_generators = 3
    power_factor = 0.9

    # Load parameters
    p_base_load_per_unit = 50.0 / 3.0  # MW per base load
    num_base_loads = 3
    p_additional_load = 100.0  # MW, active load at Bus 8

    # --- Step 1: Calculate total active power generation ---
    s_gen_total = num_generators * s_gen_per_unit
    p_gen_total = s_gen_total * power_factor

    # --- Step 2: Calculate total active power load ---
    p_base_load_total = num_base_loads * p_base_load_per_unit
    p_load_total = p_base_load_total + p_additional_load

    # --- Step 3: Calculate total power loss ---
    p_loss_total = p_gen_total - p_load_total

    # --- Step 4: Print the results and the final equation ---
    print("Calculation of Total Power Loss:")
    print("-" * 35)
    print(f"Total Active Power Generation (P_gen):")
    print(f"  {num_generators} generators * {s_gen_per_unit:.1f} MVA/gen * {power_factor} PF = {p_gen_total:.3f} MW")
    print("\nTotal Active Power Load (P_load):")
    print(f"  Base Loads: {num_base_loads} * {p_base_load_per_unit:.3f} MW = {p_base_load_total:.3f} MW")
    print(f"  Additional Load: {p_additional_load:.3f} MW")
    print(f"  Total Load = {p_base_load_total:.3f} MW + {p_additional_load:.3f} MW = {p_load_total:.3f} MW")
    print("-" * 35)
    print("\nFinal Equation for Total Power Loss:")
    print(f"P_loss = P_gen - P_load")
    print(f"P_loss = {p_gen_total:.3f} MW - {p_load_total:.3f} MW = {p_loss_total:.3f} MW")

if __name__ == '__main__':
    calculate_power_losses()