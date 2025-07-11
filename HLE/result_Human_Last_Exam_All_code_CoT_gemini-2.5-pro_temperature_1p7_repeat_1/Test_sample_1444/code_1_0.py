import math

def solve_tree_journey():
    """
    Calculates the walking tree's total distance traveled and direction.
    """
    # Data provided in the problem
    initial_angle = 45.0  # degrees, in 1000 CE
    final_angle = 34.4    # degrees, in 2000 CE
    initial_walk_m = 100  # meters

    # --- Step 1: Determine the direction ---
    # Since the angle to Polaris (latitude) decreased, the tree moved south.
    direction = "South"
    print(f"Direction of Travel: {direction}")
    print("-" * 30)

    # --- Step 2 & 3: Calculate distance from angle changes ---
    # The distance for one degree of latitude is based on Earth's polar circumference (~40,008 km).
    km_per_degree_latitude = 40008 / 360.0

    # Calculate total change in latitude
    angle_change_degrees = initial_angle - final_angle
    
    # Calculate distance from the main walks
    distance_from_angles_km = angle_change_degrees * km_per_degree_latitude
    
    print("Calculating distance from latitude change:")
    print(f"Total Latitude Change = Initial Angle ({initial_angle}°) - Final Angle ({final_angle}°)")
    print(f"                      = {angle_change_degrees:.1f}°")
    print("\nUsing the conversion factor of {:.3f} km per degree.".format(km_per_degree_latitude))
    print(f"Distance = Latitude Change ({angle_change_degrees:.1f}°) * Kilometers per Degree ({km_per_degree_latitude:.3f} km/°)")
    print(f"         = {distance_from_angles_km:.3f} km")
    print("-" * 30)

    # --- Step 4: Calculate total distance including the initial walk ---
    initial_walk_km = initial_walk_m / 1000.0
    total_distance_km = distance_from_angles_km + initial_walk_km
    
    print("Calculating total distance traveled:")
    print(f"Total Distance = Distance from Century Walks ({distance_from_angles_km:.3f} km) + Initial Walk ({initial_walk_km:.3f} km)")
    print(f"               = {total_distance_km:.3f} km")
    print("-" * 30)
    
    # --- Step 5: Format the final answer ---
    final_value = total_distance_km * 10
    final_answer = round(final_value)

    print("Final answer calculation:")
    print(f"Final Value = Total Distance in km ({total_distance_km:.3f}) * 10")
    print(f"            = {final_value:.3f}")
    print(f"Rounded to the nearest integer, the answer is: {final_answer}")
    
    return final_answer

final_answer = solve_tree_journey()
print(f"<<<{final_answer}>>>")