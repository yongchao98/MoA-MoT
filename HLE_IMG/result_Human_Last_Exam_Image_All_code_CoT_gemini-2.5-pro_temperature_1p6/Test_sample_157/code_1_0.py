import math

def calculate_power_losses():
    """
    Calculates the total resistive power losses in the power system based on the provided data.
    """

    # --- System Parameters ---
    
    # Generation
    num_generators = 3
    generator_s_rating_mva = 100.0  # Apparent power rating per generator in MVA
    system_pf = 0.9  # Power factor for all components, lagging

    # Loads
    base_load_mw = 50.0 / 3.0  # Active power of each base load in MW
    num_base_loads = 3
    additional_load_mw = 100.0  # Active power of the additional load in MW

    # --- Step 1: Calculate Total Generated Active Power (P_gen) ---
    total_generator_s_mva = num_generators * generator_s_rating_mva
    total_p_gen_mw = total_generator_s_mva * system_pf

    # --- Step 2: Calculate Total Active Power Load (P_load) ---
    total_base_load_mw = num_base_loads * base_load_mw
    total_p_load_mw = total_base_load_mw + additional_load_mw

    # --- Step 3: Calculate Total Resistive Power Loss (P_loss) ---
    total_p_loss_mw = total_p_gen_mw - total_p_load_mw
    
    # --- Output the results step-by-step as requested ---
    
    print("This script calculates the total power loss based on the principle of power balance (P_loss = P_generated - P_load).\n")
    
    print("1. Calculating Total Active Power Generated (P_gen):")
    print(f"Total Generation Capacity = {num_generators} generators * {generator_s_rating_mva} MVA/generator = {total_generator_s_mva:.3f} MVA")
    print(f"Power Factor = {system_pf}")
    print(f"P_gen = {total_generator_s_mva:.3f} MVA * {system_pf} = {total_p_gen_mw:.3f} MW\n")

    print("2. Calculating Total Active Power Load (P_load):")
    print(f"Total Base Load = {num_base_loads} loads * {base_load_mw:.3f} MW/load = {total_base_load_mw:.3f} MW")
    print(f"Additional Load = {additional_load_mw:.3f} MW")
    print(f"P_load = {total_base_load_mw:.3f} MW + {additional_load_mw:.3f} MW = {total_p_load_mw:.3f} MW\n")
    
    print("3. Calculating Total Power Loss (P_loss):")
    print("The final equation is P_loss = P_gen - P_load")
    print(f"P_loss = {total_p_gen_mw:.3f} - {total_p_load_mw:.3f} = {total_p_loss_mw:.3f} MW\n")
    
    # Final Answer
    print(f"The total power losses in the system are {total_p_loss_mw:.3f} MW.")

if __name__ == "__main__":
    calculate_power_losses()
    
    # Deriving the value for the final answer format
    total_p_loss_mw = (3 * 100.0 * 0.9) - ((50.0 / 3.0 * 3) + 100.0)
    print(f'<<<{total_p_loss_mw:.3f}>>>')