import math

def solve_walking_tree():
    """
    Calculates the total distance traveled by the Ancient Walking Tree and its direction.
    """
    # --- Data from the problem ---
    initial_lat = 45.0  # Latitude in 1000 CE
    first_walk_lat = 44.8 # Latitude in 1100 CE
    final_lat = 34.4    # Latitude in 2000 CE
    first_walk_dist_m = 100 # Distance of the first walk in meters

    # --- Determine Direction ---
    # The latitude decreases, which means the tree is moving South.
    direction = "South"
    print(f"The direction the tree was walking: {direction}")
    print("-" * 20)

    # --- Calculate Total Distance ---
    print("Calculating the approximate total distance traveled:\n")

    # Step 1: Calculate the change in latitude for the first walk
    lat_change_first_walk = initial_lat - first_walk_lat
    print("Step 1: Calculate the latitude change for the first walk to find the conversion ratio.")
    print(f"   - Change in latitude (1000-1100 CE) = {initial_lat}° - {first_walk_lat}° = {lat_change_first_walk:.1f}°")

    # Step 2: Establish the distance-to-latitude-change ratio
    dist_per_degree = first_walk_dist_m / lat_change_first_walk
    print("\nStep 2: Establish the distance-to-latitude-change ratio from the first walk.")
    print(f"   - Ratio = {first_walk_dist_m} meters / {lat_change_first_walk:.1f}° = {dist_per_degree:.1f} meters per degree")

    # Step 3: Calculate the total change in latitude over 1000 years
    total_lat_change = initial_lat - final_lat
    print("\nStep 3: Calculate the total change in latitude from 1000 CE to 2000 CE.")
    print(f"   - Total change in latitude = {initial_lat}° - {final_lat}° = {total_lat_change:.1f}°")

    # Step 4: Calculate the total distance traveled
    total_dist_m = total_lat_change * dist_per_degree
    print("\nStep 4: Calculate the total distance traveled in meters.")
    print(f"   - Total distance = {total_lat_change:.1f}° * {dist_per_degree:.1f} meters/degree = {total_dist_m:.1f} meters")

    # Step 5: Convert the total distance to kilometers
    total_dist_km = total_dist_m / 1000
    print("\nStep 5: Convert total distance to kilometers.")
    print(f"   - Total distance in km = {total_dist_m:.1f} meters / 1000 = {total_dist_km:.1f} km")
    print("-" * 20)

    # Step 6: Calculate and print the final answer in the requested format
    final_value = total_dist_km * 10
    final_answer = round(final_value)
    print("Final Answer Calculation:")
    print(f"   - The final value is Nearest Integer({total_dist_km:.1f} * 10)")
    print(f"   - Which is Nearest Integer({final_value:.1f})")
    print(f"   - Result: {final_answer}")
    
    # Final answer in the specified format
    print(f"\n<<<{final_answer}>>>")


solve_walking_tree()