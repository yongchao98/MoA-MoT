import math

def solve_tree_puzzle():
    """
    Calculates the total distance traveled by the Ancient Walking Tree and its direction.
    """
    # Data from the problem
    initial_pos_1000ce = 45.0  # degrees
    pos_after_first_walk_1100ce = 44.8  # degrees
    final_pos_2000ce = 34.4  # degrees
    distance_first_walk_m = 100  # meters

    # 1. Determine the direction
    # The latitude decreases from 45.0 to 34.4 degrees.
    # A decrease in latitude in the Northern Hemisphere means moving south.
    direction = "South"

    # 2. Establish the scale (km per degree of angle change)
    # We assume the "first walk" of 100m corresponds to the first recorded change in angle.
    distance_first_walk_km = distance_first_walk_m / 1000.0
    angle_change_first_walk = initial_pos_1000ce - pos_after_first_walk_1100ce
    
    # Avoid division by zero, though it's not an issue here
    if angle_change_first_walk == 0:
        print("Error: Angle change for the first walk is zero, cannot determine scale.")
        return

    scaling_factor_km_per_degree = distance_first_walk_km / angle_change_first_walk

    # 3. Calculate the total distance traveled
    total_angle_change = initial_pos_1000ce - final_pos_2000ce
    total_distance_km = scaling_factor_km_per_degree * total_angle_change

    # 4. Calculate the final answer as per the required format
    final_answer_value = total_distance_km * 10
    final_answer = round(final_answer_value)

    # --- Output the results ---
    print("Step 1: Determine the Direction of Travel")
    print(f"The tree's latitude changed from {initial_pos_1000ce}° to {final_pos_2000ce}°. This is a southward movement.")
    print(f"Direction: {direction}\n")

    print("Step 2: Calculate the Scaling Factor")
    print(f"The first walk resulted in an angle change of {initial_pos_1000ce}° - {pos_after_first_walk_1100ce}° = {angle_change_first_walk:.1f}°.")
    print(f"This walk covered {distance_first_walk_m} meters, or {distance_first_walk_km} km.")
    print(f"Scaling Factor = {distance_first_walk_km} km / {angle_change_first_walk:.1f}° = {scaling_factor_km_per_degree:.1f} km/degree\n")

    print("Step 3: Calculate the Total Distance")
    print(f"Total Angle Change = {initial_pos_1000ce}° - {final_pos_2000ce}° = {total_angle_change:.1f}°")
    print("Final Equation:")
    print(f"Total Distance (km) = Scaling Factor * Total Angle Change")
    print(f"Total Distance (km) = {scaling_factor_km_per_degree:.1f} * {total_angle_change:.1f} = {total_distance_km:.1f} km\n")
    
    print("The approximate total distance the tree has traveled is 5.3 km.")
    print("The direction the tree was walking is South.")

    print(f"\nFinal Answer Calculation: round({total_distance_km:.1f} * 10) = round({final_answer_value:.1f}) = {final_answer}")
    print(f"<<<{final_answer}>>>")

solve_tree_puzzle()