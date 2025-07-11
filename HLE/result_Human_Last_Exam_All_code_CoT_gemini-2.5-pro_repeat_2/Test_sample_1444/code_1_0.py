import math

def solve_walking_tree():
    """
    Calculates the total distance and direction of the Ancient Walking Tree's travel.
    """
    # Data provided in the problem
    positions = {
        1000: 45.0, 1100: 44.8, 1200: 44.3, 1300: 43.5, 1400: 42.4,
        1500: 41.0, 1600: 39.4, 1700: 37.8, 1800: 36.5, 1900: 35.2,
        2000: 34.4
    }
    distance_first_walk_m = 100

    # Step 1: Determine the direction
    print("Step 1: Determining the Direction of Travel")
    initial_pos = positions[1000]
    final_pos = positions[2000]
    direction = "South" if final_pos < initial_pos else "North"
    print(f"The angle towards Polaris decreased from {initial_pos}° in 1000 CE to {final_pos}° in 2000 CE.")
    print(f"A decrease in the angle to the North Star corresponds to a decrease in latitude.")
    print(f"Therefore, the direction the tree was walking is: {direction}\n")

    # Step 2: Calculate the distance-to-angle scaling factor from the first walk
    print("Step 2: Calculating the Scaling Factor")
    pos_1000 = positions[1000]
    pos_1100 = positions[1100]
    angular_change_first_walk = pos_1000 - pos_1100
    distance_first_walk_km = distance_first_walk_m / 1000.0
    scaling_factor_km_per_degree = distance_first_walk_km / angular_change_first_walk
    
    print(f"The angular change during the first walk (1000 CE to 1100 CE) was: {pos_1000}° - {pos_1100}° = {angular_change_first_walk:.1f}°")
    print(f"The distance covered in this walk was {distance_first_walk_m} m, which is {distance_first_walk_km} km.")
    print(f"This establishes a scaling factor: {distance_first_walk_km} km / {angular_change_first_walk:.1f}° = {scaling_factor_km_per_degree} km per degree.\n")

    # Step 3: Calculate the total distance traveled
    print("Step 3: Calculating the Total Distance Traveled")
    total_angular_change = initial_pos - final_pos
    total_distance_km = total_angular_change * scaling_factor_km_per_degree
    
    print("The total distance is calculated by the formula:")
    print("Total Distance = (Total Angular Change) * (Scaling Factor)")
    print(f"The numbers for the final equation are:")
    print(f"Total Distance = ({initial_pos} - {final_pos}) * {scaling_factor_km_per_degree}")
    print(f"Total Distance = {total_angular_change:.1f}° * {scaling_factor_km_per_degree} km/°")
    print(f"The approximate total distance the tree has traveled is: {total_distance_km:.1f} km\n")

    # Step 4: Calculate the final answer in the required format
    final_answer = round(total_distance_km * 10)
    print("The required final answer is: Nearest Integer(Total Distance in km * 10)")
    print(f"Calculation: round({total_distance_km:.1f} * 10) = round({total_distance_km * 10}) = {final_answer}")
    
    # Final answer in the specified format
    print(f"\n<<<53>>>")

solve_walking_tree()