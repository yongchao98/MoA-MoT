import math

def solve_walking_tree_problem():
    """
    Calculates the total distance and direction traveled by the Ancient Walking Tree
    and provides the final answer in the required format.
    """

    # Data from the problem
    positions = {
        1000: 45.0, 1100: 44.8, 1200: 44.3, 1300: 43.5, 1400: 42.4,
        1500: 41.0, 1600: 39.4, 1700: 37.8, 1800: 36.5, 1900: 35.2,
        2000: 34.4
    }
    first_walk_dist_m = 100.0

    # --- Part 1: Determine the Direction ---
    print("Part 1: Direction of Travel")
    print("==============================")
    print("The angle of Polaris above the horizon is equal to the observer's latitude.")
    print(f"The tree's latitude decreased from {positions[1000]}° in 1000 CE to {positions[2000]}° in 2000 CE.")
    print("In the Northern Hemisphere, a decrease in latitude means moving south.")
    print("Therefore, the tree was walking directly South.\n")

    # --- Part 2: Calculate the Total Distance ---
    print("Part 2: Total Distance Traveled")
    print("=================================")
    print("To find the total distance, we first need to establish a scale from the given data.")

    # Step 1: Calibration
    lat_initial = positions[1000]
    lat_after_first_walk = positions[1100]
    lat_final = positions[2000]

    first_walk_lat_change = lat_initial - lat_after_first_walk
    first_walk_dist_km = first_walk_dist_m / 1000.0
    
    print(f"The first walk resulted in a latitude change of {lat_initial}° - {lat_after_first_walk}° = {first_walk_lat_change:.1f}°.")
    print(f"This change in latitude corresponds to a known distance of {int(first_walk_dist_m)} meters.")

    # Step 2: Calculate conversion factor (km per degree)
    km_per_degree = first_walk_dist_km / first_walk_lat_change
    print(f"This gives us a conversion factor: {first_walk_dist_km} km / {first_walk_lat_change:.1f}° = {km_per_degree:.1f} km per degree.\n")

    # Step 3: Calculate total distance
    total_lat_change = lat_initial - lat_final
    total_distance_km = total_lat_change * km_per_degree

    print("Now, we can calculate the total distance traveled over 1000 years.")
    print(f"The total latitude change is {lat_initial}° - {lat_final}° = {total_lat_change:.1f}°.")
    
    print("\nThe final equation for the total distance in km is:")
    print(f"Total Distance (km) = (Total Latitude Change) * (km per degree)")
    print(f"Total Distance (km) = {total_lat_change:.1f} * {km_per_degree:.1f} = {total_distance_km:.1f} km\n")

    # --- Final Answer Formatting ---
    print("Final Answer Calculation")
    print("========================")
    final_value_float = total_distance_km * 10
    final_answer = round(final_value_float)
    
    print(f"The required format is: Nearest Integer(Total Distance in km * 10)")
    print(f"Calculation: round({total_distance_km:.1f} * 10) = round({final_value_float:.1f}) = {final_answer}")
    
    # Output the final answer in the required format
    print(f"<<<{final_answer}>>>")

solve_walking_tree_problem()