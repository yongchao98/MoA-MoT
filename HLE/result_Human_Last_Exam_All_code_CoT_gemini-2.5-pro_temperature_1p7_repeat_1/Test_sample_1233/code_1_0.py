import math

def simulate_growth():
    """
    This simulation demonstrates how extending incubation time can help in identifying
    a slow-growing pathogen like Campylobacter when co-cultured with a
    fast-growing organism like Bacillus.
    """
    
    # --- Parameters ---
    # Assume a very small initial colony size (e.g., a single cell)
    initial_colony_size = 0.01  # in mm
    
    # Growth rates (doubling time in hours)
    # Bacillus species are often fast growers.
    bacillus_doubling_time = 6  # hours to double in size
    # Campylobacter is a slow grower.
    campylobacter_doubling_time = 18 # hours to double in size
    
    # Visibility threshold in mm
    visibility_threshold = 0.5 # mm
    
    # Incubation times to check
    first_check_time = 48  # hours (2 days)
    second_check_time = 72 # hours (3 days)
    
    print("--- Lab Simulation: Bacterial Growth on Plate ---\n")
    print(f"Investigating plates at the {first_check_time}-hour mark (Day 2):\n")
    
    # --- Calculations for 48 hours ---
    # Growth formula: final_size = initial_size * 2^(time / doubling_time)
    num_doubles_bacillus_t1 = first_check_time / bacillus_doubling_time
    size_bacillus_t1 = initial_colony_size * (2 ** num_doubles_bacillus_t1)
    
    num_doubles_campy_t1 = first_check_time / campylobacter_doubling_time
    size_campy_t1 = initial_colony_size * (2 ** num_doubles_campy_t1)
    
    print(f"Bacillus colony size: {initial_colony_size} * 2^({first_check_time}/{bacillus_doubling_time}) = {size_bacillus_t1:.2f} mm")
    print(f"Campylobacter colony size: {initial_colony_size} * 2^({first_check_time}/{campylobacter_doubling_time}) = {size_campy_t1:.2f} mm")
    
    if size_bacillus_t1 > visibility_threshold:
        print("Result: Large Bacillus colonies are clearly visible.")
    if size_campy_t1 < visibility_threshold:
        print("Result: Campylobacter colonies are too small to be seen and could be missed.")
        
    print("\n-------------------------------------------------\n")
    print(f"Action: Incubate for longer. Re-investigating at the {second_check_time}-hour mark (Day 3):\n")
    
    # --- Calculations for 72 hours ---
    num_doubles_bacillus_t2 = second_check_time / bacillus_doubling_time
    size_bacillus_t2 = initial_colony_size * (2 ** num_doubles_bacillus_t2)
    
    num_doubles_campy_t2 = second_check_time / campylobacter_doubling_time
    size_campy_t2 = initial_colony_size * (2 ** num_doubles_campy_t2)
    
    print(f"Bacillus colony size: {initial_colony_size} * 2^({second_check_time}/{bacillus_doubling_time}) = {size_bacillus_t2:.2f} mm (likely overgrown)")
    print(f"Final Campylobacter size calculation: {initial_colony_size} * 2^({second_check_time} / {campylobacter_doubling_time}) = {size_campy_t2:.2f} mm")

    if size_campy_t2 > visibility_threshold:
        print("Result: Campylobacter colonies have now grown to a visible size and can be identified.")

simulate_growth()