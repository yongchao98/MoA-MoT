import math

def solve_ancient_tree_problem():
    """
    Calculates the total distance and direction traveled by the Ancient Walking Tree.
    """
    initial_latitude = 45.0  # degrees, in 1000 CE before the first walk
    final_latitude = 34.4    # degrees, in 2000 CE after the last walk
    
    # One degree of latitude is approximately 111.32 km.
    # This is based on Earth's circumference of ~40,075 km (40075 / 360 = 111.319...).
    km_per_degree_latitude = 111.32

    # --- Step 1: Determine Direction ---
    # The latitude decreases, so the tree is moving South.
    direction = "South"
    
    # --- Step 2: Calculate Total Distance ---
    # Calculate the total change in latitude
    latitude_change = initial_latitude - final_latitude
    
    # Calculate the total distance in kilometers
    total_distance_km = latitude_change * km_per_degree_latitude

    # --- Step 3: Present the Results and Final Answer ---
    print("Direction of Travel:")
    print(f"The tree's latitude decreased from {initial_latitude}° to {final_latitude}°, so it walked {direction}.")
    print("\nTotal Distance Calculation:")
    print(f"Total latitude change = {initial_latitude}° - {final_latitude}° = {latitude_change:.1f}°")
    print(f"Using the approximation of {km_per_degree_latitude} km per degree of latitude:")
    print(f"Total distance (km) = {latitude_change:.1f} * {km_per_degree_latitude} = {total_distance_km:.3f} km")
    
    # The problem asks for the nearest integer of (Total Distance in km * 10)
    final_value_raw = total_distance_km * 10
    final_answer = round(final_value_raw)
    
    print("\nFinal Answer Formulation:")
    print(f"The required value is Nearest Integer(Total Distance in km * 10).")
    print(f"Calculation: round({total_distance_km:.3f} * 10) = round({final_value_raw:.3f}) = {final_answer}")

    # Final answer in the required format
    print(f"\n<<<>>>")
    print(f"<<<{final_answer}>>>")
    print(f"<<<>>>")

solve_ancient_tree_problem()