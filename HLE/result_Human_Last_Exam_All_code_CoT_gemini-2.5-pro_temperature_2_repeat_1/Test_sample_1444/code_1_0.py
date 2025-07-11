import math

def solve_ancient_tree_mystery():
    """
    Calculates the total distance traveled by the Ancient Walking Tree and its direction.
    """
    # Recorded positions (latitude in degrees) of the tree each century
    positions = {
        1000: 45.0, 1100: 44.8, 1200: 44.3, 1300: 43.5,
        1400: 42.4, 1500: 41.0, 1600: 39.4, 1700: 37.8,
        1800: 36.5, 1900: 35.2, 2000: 34.4
    }
    
    initial_lat = positions[1000]
    final_lat = positions[2000]

    # --- Step 1: Determine the Direction ---
    print("Step 1: Determining the Direction of Travel")
    print(f"The tree's latitude changed from {initial_lat}° in 1000 CE to {final_lat}° in 2000 CE.")
    print("Since the latitude (angle to Polaris) decreased, the tree has been consistently walking South.")
    direction = "South"
    print(f"Direction of travel: {direction}\n")
    
    # --- Step 2: Establish the Distance-to-Angle Scale ---
    print("Step 2: Establishing the Distance Scale")
    lat_1000 = positions[1000]
    lat_1100 = positions[1100]
    distance_first_walk_m = 100  # in meters

    # The change in latitude during the first walk
    lat_change_first_walk = lat_1000 - lat_1100
    
    # The scale in meters per degree
    m_per_degree = distance_first_walk_m / lat_change_first_walk
    
    print(f"The first walk was {distance_first_walk_m} meters long.")
    print(f"The latitude changed by {lat_1000}° - {lat_1100}° = {lat_change_first_walk:.1f}°.")
    print(f"This gives us a scale of {distance_first_walk_m} m / {lat_change_first_walk:.1f}° = {m_per_degree:.0f} meters per degree.\n")
    
    # --- Step 3: Calculate the Total Distance ---
    print("Step 3: Calculating Total Distance")
    # Total change in latitude over 1000 years
    total_lat_change = initial_lat - final_lat
    
    # Total distance in meters
    total_distance_m = total_lat_change * m_per_degree
    total_distance_km = total_distance_m / 1000

    print(f"Total latitude change = {initial_lat}° - {final_lat}° = {total_lat_change:.1f}°.")
    print("Final Equation for Total Distance in meters:")
    print(f"Total Distance (m) = Total Latitude Change * Meters Per Degree")
    print(f"Total Distance (m) = {total_lat_change:.1f} * {m_per_degree:.0f} = {total_distance_m:.0f} m")
    print(f"This is equal to {total_distance_km:.1f} km.\n")
    
    # --- Step 4: Final Answer Calculation ---
    print("Step 4: Final Answer Calculation")
    final_value = total_distance_km * 10
    final_answer = round(final_value)
    
    print("Final calculation: Nearest Integer(Total Distance in km * 10)")
    print(f"Value = Nearest Integer({total_distance_km:.1f} * 10)")
    print(f"Value = Nearest Integer({final_value:.1f})")
    print(f"The final answer is {final_answer}.")
    
    # --- Final Result ---
    print("\n--- Summary ---")
    print(f"Approximate total distance traveled: {total_distance_km:.1f} km")
    print(f"Direction of walking: {direction}")
    
    # Adhering to the specified output format
    print(f"\n<<<{final_answer}>>>")

solve_ancient_tree_mystery()