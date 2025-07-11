import math

def solve_tree_problem():
    """
    Calculates the total distance traveled by the walking tree and its direction.
    """
    # --- Data from the problem ---
    initial_lat = 45.0
    lat_after_first_walk = 44.8
    final_lat = 34.4
    dist_first_walk_meters = 100.0

    # --- 1. Determine the direction ---
    # In the Northern Hemisphere, latitude equals the angle of Polaris.
    # A decrease in latitude means moving South.
    direction = "South"
    print(f"Direction of Travel: {direction}")
    print("-" * 40)

    # --- 2. Calculate the distance scale ---
    lat_change_first_walk = initial_lat - lat_after_first_walk
    meters_per_degree = dist_first_walk_meters / lat_change_first_walk

    print("Calculating the distance-to-latitude scale:")
    print(f"The first walk changed latitude by: {initial_lat}° - {lat_after_first_walk}° = {lat_change_first_walk:.1f}°")
    print(f"This covered {dist_first_walk_meters:.0f} meters.")
    print(f"Scale = {dist_first_walk_meters:.0f} meters / {lat_change_first_walk:.1f}° = {meters_per_degree:.0f} meters per degree")
    print("-" * 40)

    # --- 3. Calculate total distance ---
    total_lat_change = initial_lat - final_lat
    total_dist_meters = total_lat_change * meters_per_degree
    total_dist_km = total_dist_meters / 1000.0

    print("Calculating the total distance traveled from 1000 CE to 2000 CE:")
    print("The equation for the total distance is:")
    print(f"Distance (m) = (Initial Latitude - Final Latitude) * Scale")
    print(f"Distance (m) = ({initial_lat} - {final_lat}) * {meters_per_degree:.0f}")
    print(f"Distance (m) = {total_lat_change:.1f} * {meters_per_degree:.0f} = {total_dist_meters:.0f} meters")
    print(f"Total distance in kilometers: {total_dist_km:.1f} km")
    print("-" * 40)

    # --- 4. Final Answer Formatting ---
    final_answer = round(total_dist_km * 10)
    print("Final answer calculation:")
    print(f"Result = Nearest Integer(Total Distance in km * 10)")
    print(f"Result = round({total_dist_km:.1f} * 10) = round({total_dist_km * 10:.1f}) = {final_answer}")
    
    # --- Required Output Format ---
    print(f"\n<<<{final_answer}>>>")

solve_tree_problem()