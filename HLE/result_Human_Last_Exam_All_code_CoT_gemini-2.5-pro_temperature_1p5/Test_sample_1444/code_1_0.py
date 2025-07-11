import math

def solve_walking_tree():
    """
    Solves the Ancient Walking Tree problem.
    """
    # Latitude readings from 1000 CE to 2000 CE
    latitudes = [45.0, 44.8, 44.3, 43.5, 42.4, 41.0, 39.4, 37.8, 36.5, 35.2, 34.4]
    
    # 1. Establish the conversion factor from the first walk
    # The change in latitude for the first walk was 45.0 - 44.8 = 0.2 degrees
    # This walk covered 100 meters.
    first_walk_lat_change = latitudes[0] - latitudes[1]
    first_walk_distance_m = 100.0
    conversion_factor_m_per_deg = first_walk_distance_m / first_walk_lat_change
    
    # 2. Calculate the distance of each walk and the total distance
    total_distance_m = 0
    individual_distances_m = []
    
    for i in range(len(latitudes) - 1):
        lat_start = latitudes[i]
        lat_end = latitudes[i+1]
        
        # The change in latitude is always positive in our calculation (absolute change)
        lat_change = abs(lat_end - lat_start)
        
        # Calculate distance for this specific walk
        walk_distance_m = lat_change * conversion_factor_m_per_deg
        individual_distances_m.append(walk_distance_m)
        total_distance_m += walk_distance_m

    # 3. Print the results and the breakdown of the calculation
    print("Calculating the total distance traveled by the tree:")
    
    # Format the equation string
    # Using round() to handle potential floating point inaccuracies for cleaner display
    equation_parts = [str(int(round(d))) for d in individual_distances_m]
    equation_str = " + ".join(equation_parts)
    
    print(f"Total Distance (meters) = {equation_str}")
    print(f"Total Distance (meters) = {int(round(total_distance_m))} meters")
    
    total_distance_km = total_distance_m / 1000.0
    print(f"Total Distance (km) = {total_distance_km:.1f} km")
    
    # 4. Determine and print the direction
    direction = "South"
    print(f"\nDirection of travel: {direction}")
    print("This is because the latitude (angle to Polaris) consistently decreases.")
    
    # 5. Calculate and format the final answer
    final_answer_val = round(total_distance_km * 10)
    
    print("\nFinal calculation as per the question:")
    print(f"Nearest Integer(Total Distance in km * 10) = round({total_distance_km:.1f} * 10)")
    print(f"Result = round({total_distance_km * 10:.1f}) = {final_answer_val}")
    
    # Final answer in the specified format
    print(f"\n<<< {final_answer_val} >>>")

solve_walking_tree()