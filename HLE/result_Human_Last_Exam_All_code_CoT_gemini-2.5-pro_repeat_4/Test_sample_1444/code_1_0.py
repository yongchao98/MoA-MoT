import math

def solve_ancient_tree_problem():
    """
    Calculates the total distance and direction traveled by the Ancient Walking Tree.
    """
    # Position data: Year -> Angle of Polaris (in degrees)
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
        2000: 34.4
    }

    # --- Part 1: Determine the direction ---
    initial_angle = positions[1000]
    final_angle = positions[2000]

    if final_angle < initial_angle:
        direction = "South"
    elif final_angle > initial_angle:
        direction = "North"
    else:
        direction = "Stationary (East/West or no movement)"

    print(f"The direction the tree was walking: {direction}\n")

    # --- Part 2: Calculate the total distance ---
    print("Calculating the total distance traveled:")

    # Calibration using the first walk
    first_walk_dist_m = 100.0
    first_walk_angle_change = positions[1000] - positions[1100]
    
    # This is the conversion factor specific to this problem
    meters_per_degree = first_walk_dist_m / first_walk_angle_change

    # Calculate total change in angle
    total_angle_change = initial_angle - final_angle
    
    # Calculate total distance
    total_distance_m = total_angle_change * meters_per_degree
    total_distance_km = total_distance_m / 1000

    print("Step 1: Determine the conversion factor from the first walk.")
    print(f"  - Angle change: {positions[1000]:.1f}° - {positions[1100]:.1f}° = {first_walk_angle_change:.1f}°")
    print(f"  - Distance covered: {int(first_walk_dist_m)} meters")
    print(f"  - Factor: {int(first_walk_dist_m)} meters / {first_walk_angle_change:.1f}° = {int(meters_per_degree)} meters per degree\n")

    print("Step 2: Calculate the total distance over 1000 years.")
    print(f"  - Total angle change: {initial_angle:.1f}° - {final_angle:.1f}° = {total_angle_change:.1f}°")
    print(f"  - Total distance (m): {total_angle_change:.1f}° * {int(meters_per_degree)} m/° = {total_distance_m:.0f} meters")
    print(f"  - Total distance (km): {total_distance_km:.1f} km\n")
    
    # --- Part 3: Format the final answer ---
    final_answer = round(total_distance_km * 10)
    
    print("Final answer format: Nearest Integer(Total Distance in km * 10)")
    print(f"Calculation: round({total_distance_km:.1f} * 10) = round({total_distance_km * 10:.1f}) = {final_answer}")
    
    print(f"<<<{final_answer}>>>")

solve_ancient_tree_problem()