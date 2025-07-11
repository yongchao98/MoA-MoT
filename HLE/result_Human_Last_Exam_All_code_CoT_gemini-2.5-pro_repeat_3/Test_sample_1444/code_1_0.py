import math

def solve_ancient_tree_problem():
    """
    Calculates the total distance traveled and direction of the Ancient Walking Tree.
    """
    # Position data: Year -> Latitude (degrees)
    positions = {
        1000: 45.0,
        1100: 44.8,
        1200: 44.3,
        1300: 43.5,
        1400: 42.4,
        1500: 41.0,
        1600: 39.4,
        1700: 37.8,
        1800: 36.5,
        1900: 35.2,
        2000: 34.4,
    }

    # Part 1: Determine the direction
    # The angle of Polaris is the latitude. Since it's decreasing, the tree is moving south.
    initial_lat = positions[1000]
    final_lat = positions[2000]
    direction = "South" if final_lat < initial_lat else "North"
    
    print(f"The tree's latitude decreased from {initial_lat}° to {final_lat}°. Therefore, the direction of walking is {direction}.")
    print("-" * 20)

    # Part 2: Calculate the total distance
    
    # Step 2a: Find the conversion factor from the first walk (1000 CE to 1100 CE)
    lat_change_first_walk = positions[1000] - positions[1100]
    dist_first_walk_m = 100
    dist_first_walk_km = dist_first_walk_m / 1000.0
    
    # km per degree of latitude change
    km_per_degree = dist_first_walk_km / lat_change_first_walk
    
    print("Calculating the distance conversion factor:")
    print(f"The first walk resulted in a latitude change of {positions[1000]:.1f}° - {positions[1100]:.1f}° = {lat_change_first_walk:.1f}°.")
    print(f"This change corresponds to a distance of {dist_first_walk_m} m, or {dist_first_walk_km} km.")
    print(f"Conversion Factor = {dist_first_walk_km} km / {lat_change_first_walk:.1f}° = {km_per_degree:.1f} km per degree.")
    print("-" * 20)
    
    # Step 2b: Calculate the total distance traveled from 1000 CE to 2000 CE
    total_lat_change = initial_lat - final_lat
    total_distance_km = total_lat_change * km_per_degree
    
    print("Calculating the total distance traveled:")
    print(f"Total latitude change = {initial_lat:.1f}° - {final_lat:.1f}° = {total_lat_change:.1f}°.")
    print("The final equation for the total distance is:")
    print(f"Total Distance (km) = {total_lat_change:.1f} * {km_per_degree:.1f} = {total_distance_km:.1f} km")
    print("-" * 20)
    
    # Part 3: Calculate and format the final answer
    final_value = total_distance_km * 10
    final_answer = round(final_value)
    
    print("Final answer calculation:")
    print(f"Value = Total Distance (km) * 10 = {total_distance_km:.1f} * 10 = {final_value:.1f}")
    print(f"Rounding to the nearest integer gives: {final_answer}")
    
    print(f"\n<<< {final_answer} >>>")

solve_ancient_tree_problem()