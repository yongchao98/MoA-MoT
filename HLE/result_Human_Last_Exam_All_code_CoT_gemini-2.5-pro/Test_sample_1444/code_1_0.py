import math

def solve_ancient_tree_problem():
    """
    Solves the Ancient Walking Tree problem by calculating the total distance
    and direction of its travel over 1000 years.
    """
    # Data from the problem description
    positions = {
        1000: 45.0, 1100: 44.8, 1200: 44.3, 1300: 43.5, 1400: 42.4,
        1500: 41.0, 1600: 39.4, 1700: 37.8, 1800: 36.5, 1900: 35.2,
        2000: 34.4
    }
    dist_first_walk_m = 100

    # Part 1: Determine the direction of travel
    # The latitude (angle relative to Polaris) decreases from 45.0 to 34.4 degrees.
    # In the Northern Hemisphere, a decrease in latitude means traveling South.
    direction = "South"
    print(f"Direction the tree was walking: {direction}")

    # Part 2: Calculate the approximate total distance traveled
    print("\nCalculating the total distance:")

    # Step 1: Establish the relationship (conversion factor) between latitude change and distance
    # using the data from the first walk (from 1000 CE to 1100 CE).
    lat_change_first_walk = positions[1000] - positions[1100]
    dist_first_walk_km = dist_first_walk_m / 1000.0
    
    # This is the km-per-degree conversion factor for this world
    km_per_degree = dist_first_walk_km / lat_change_first_walk
    
    print("\nStep 1: Find the conversion factor from latitude change to distance.")
    print(f"The first walk resulted in a latitude change of {positions[1000]:.1f}° - {positions[1100]:.1f}° = {lat_change_first_walk:.1f}°")
    print(f"This walk covered a known distance of {dist_first_walk_m} m, which is {dist_first_walk_km:.1f} km.")
    print(f"This gives a conversion factor: {dist_first_walk_km:.1f} km / {lat_change_first_walk:.1f}° = {km_per_degree:.1f} km per degree.")

    # Step 2: Calculate the total latitude change over the 1000 years.
    total_lat_change = positions[1000] - positions[2000]
    
    print("\nStep 2: Find the total change in latitude from 1000 CE to 2000 CE.")
    print(f"Total latitude change = Initial Position - Final Position")
    print(f"Total latitude change = {positions[1000]:.1f}° - {positions[2000]:.1f}° = {total_lat_change:.1f}°")

    # Step 3: Calculate the total distance traveled using the conversion factor.
    total_distance_km = total_lat_change * km_per_degree
    
    print("\nStep 3: Calculate the total distance traveled.")
    print(f"Total Distance = Total Latitude Change × Conversion Factor")
    print(f"Total Distance = {total_lat_change:.1f}° × {km_per_degree:.1f} km/° = {total_distance_km:.1f} km")
    
    print(f"\nThe approximate total distance the tree has traveled is {total_distance_km:.1f} km.")

    # Final answer calculation as requested by the format
    final_answer = round(total_distance_km * 10)
    
    print(f"<<<{final_answer}>>>")

solve_ancient_tree_problem()