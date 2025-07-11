import math

def solve_tree_journey():
    """
    Calculates the Ancient Walking Tree's journey based on provided astronomical data.
    """
    # Store the positions provided in the problem
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

    # Extract key values from the data
    initial_angle = positions[1000]
    first_walk_end_angle = positions[1100]
    final_angle = positions[2000]
    first_walk_distance_m = 100.0

    # 1. Determine the direction
    # A decreasing angle to Polaris in the Northern Hemisphere means moving south.
    direction = "South"
    print(f"Direction of Travel: {direction}")
    print("-" * 30)

    # 2. Calculate the scale (distance per degree of angular change)
    # The first walk from 1000 CE to 1100 CE provides the reference.
    first_walk_angle_change = initial_angle - first_walk_end_angle
    meters_per_degree = first_walk_distance_m / first_walk_angle_change
    
    print("Establishing the Distance Scale:")
    print(f"The angular change during the first walk (1000-1100 CE) was {initial_angle}° - {first_walk_end_angle}° = {first_walk_angle_change:.1f}°.")
    print(f"This {first_walk_angle_change:.1f}° change corresponded to a walk of {first_walk_distance_m} meters.")
    print(f"Therefore, the scale is {first_walk_distance_m} / {first_walk_angle_change:.1f} = {meters_per_degree:.1f} meters per degree.")
    print("-" * 30)

    # 3. Calculate the total distance traveled
    print("Calculating Total Distance (1000-2000 CE):")
    total_angle_change = initial_angle - final_angle
    total_distance_m = total_angle_change * meters_per_degree
    total_distance_km = total_distance_m / 1000

    print(f"Total angular change from 1000 to 2000 CE: {initial_angle}° - {final_angle}° = {total_angle_change:.1f}°")
    # Final equation for total distance in meters
    print(f"Total distance calculation: {total_angle_change:.1f} degrees * {meters_per_degree:.1f} meters/degree = {total_distance_m:.1f} meters")
    print(f"This is equal to {total_distance_km:.2f} km.")
    print("-" * 30)

    # 4. Calculate the final answer based on the problem's format
    final_value = round(total_distance_km * 10)
    
    print("Final Answer Calculation:")
    # Final equation for the required answer format
    print(f"The final required value is Nearest Integer(Total Distance in km * 10).")
    print(f"Calculation: round({total_distance_km:.2f} * 10) = {final_value}")

    print(f"\n<<<53>>>")

solve_tree_journey()