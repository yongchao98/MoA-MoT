def solve_tree_puzzle():
    """
    Solves the puzzle of the Ancient Walking Tree by calculating the total distance
    and direction of its travel, and formatting the final answer as requested.
    """
    # Store the position data provided in the problem
    positions = {
        1000: 45.0,
        1100: 44.8,
        2000: 34.4
    }
    first_walk_distance_m = 100

    # Extract key values for calculation
    initial_pos = positions[1000]
    first_walk_end_pos = positions[1100]
    final_pos = positions[2000]

    # --- Step 1: Determine the direction of travel ---
    # The angle (latitude) is decreasing, so the tree is moving South.
    direction = "South"
    print(f"1. Direction of Travel:")
    print(f"The tree's position relative to Polaris decreased from {initial_pos}° to {final_pos}°. This means its latitude decreased, so the tree walked {direction}.")
    print("-" * 30)

    # --- Step 2 & 3: Calculate the scale and total distance ---
    print("2. Total Distance Calculation:")
    # Calculate the angular change during the first walk
    first_walk_angular_change = initial_pos - first_walk_end_pos

    # Calculate the scale factor in meters per degree
    scale_factor_m_per_degree = first_walk_distance_m / first_walk_angular_change

    print(f"First, we establish a scale from the first walk (1000-1100 CE):")
    print(f"A change of ({initial_pos}° - {first_walk_end_pos}°) = {first_walk_angular_change:.1f}° corresponded to {first_walk_distance_m} meters.")
    print(f"This gives a scale factor of {scale_factor_m_per_degree:.0f} meters per degree.\n")
    
    # Calculate the total angular change over 1000 years
    total_angular_change = initial_pos - final_pos

    # Calculate the total distance in meters using the scale factor
    total_distance_m = total_angular_change * scale_factor_m_per_degree
    
    # Convert total distance to kilometers
    total_distance_km = total_distance_m / 1000

    print("Next, we calculate the total distance using the full angular change:")
    print("The final equation for the total distance is:")
    print(f"Total Distance (m) = (Initial Angle - Final Angle) * (Distance of First Walk / Angular Change of First Walk)")
    print(f"Total Distance (m) = ({initial_pos} - {final_pos}) * ({first_walk_distance_m} / ({initial_pos} - {first_walk_end_pos}))")
    print(f"Total Distance = {total_angular_change:.1f} degrees * {scale_factor_m_per_degree:.0f} m/degree = {total_distance_m:.0f} meters.")
    print(f"This is equal to {total_distance_km:.1f} kilometers.")
    print("-" * 30)

    # --- Step 4: Format the final answer ---
    final_answer = round(total_distance_km * 10)
    
    print("3. Final Answer Formatting:")
    print("The final answer is calculated as: Nearest Integer(Total Distance in km * 10)")
    print(f"Calculation: Nearest Integer({total_distance_km:.1f} * 10) = Nearest Integer({total_distance_km * 10}) = {final_answer}")

# Execute the function to print the solution
solve_tree_puzzle()
<<<53>>>