def solve_horse_problem(n, m):
    """
    Calculates the maximum water left after a horse travels m km.

    Args:
        n (int): The number of 100-liter loads of water at the origin.
        m (float): The distance to the destination in kilometers.
    """
    initial_water = n * 100.0

    # The problem assumes n*100 > m, so this check is for robustness.
    if initial_water < m:
        print("Error: Initial water is less than the distance.")
        print("It's impossible to reach the destination.")
        return

    print(f"Initial water: {initial_water:.2f} liters")
    print(f"Destination distance: {m:.2f} km")
    print("-" * 40)
    print("Strategy: Create depots to reduce trips for subsequent travel.")
    print("-" * 40)

    # --- Calculation ---
    dist_covered = 0.0
    stages_info = []  # Stores details for each stage of the journey

    # Loop for multi-load stages (from n down to 2 loads)
    for loads in range(n, 1, -1):
        if dist_covered >= m:
            break

        # For 'loads' number of 100L units, we need 2*loads - 1 trips
        trips = 2 * loads - 1
        
        # Max distance for a full stage is by using one 100L depot
        stage_dist_full = 100.0 / trips
        dist_to_destination = m - dist_covered

        if stage_dist_full <= dist_to_destination:
            # Complete this full stage
            dist_this_stage = stage_dist_full
            water_this_stage = 100.0  # By definition, a full stage consumes 100L
            stages_info.append({
                'loads': loads,
                'dist': dist_this_stage,
                'water': water_this_stage,
                'trips': trips
            })
            dist_covered += dist_this_stage
        else:
            # This is the final, partial stage before the last load
            dist_this_stage = dist_to_destination
            water_this_stage = dist_this_stage * trips
            stages_info.append({
                'loads': loads,
                'dist': dist_this_stage,
                'water': water_this_stage,
                'trips': trips
            })
            dist_covered += dist_this_stage
            # Destination is reached, no more stages needed

    # Final single-load leg, if destination not yet reached
    if dist_covered < m:
        dist_this_stage = m - dist_covered
        # With one load, there is 1 trip, so water consumed equals distance
        water_this_stage = dist_this_stage
        stages_info.append({
            'loads': 1,
            'dist': dist_this_stage,
            'water': water_this_stage,
            'trips': 1
        })

    # --- Output the results ---
    print("Calculation Breakdown of Water Consumed:")
    total_water_consumed = 0
    equation_parts = []
    
    for i, stage in enumerate(stages_info):
        print(f"\nStage {i+1} (Moving {stage['loads']} loads):")
        print(f"  Distance Covered in this Stage: {stage['dist']:.2f} km")
        print(f"  Number of Trips across this segment: {stage['trips']}")
        # The following line explicitly shows each number in the equation for this stage
        print(f"  Water Consumed in this Stage = {stage['dist']:.2f} km * {stage['trips']} trips = {stage['water']:.2f} liters")
        total_water_consumed += stage['water']
        equation_parts.append(f"{stage['water']:.2f}")

    print("-" * 40)
    final_water_left = initial_water - total_water_consumed

    # Print the final equation as requested
    print("Final Equation Summary:")
    print(f"Maximum Water Left = Initial Water - Total Water Consumed")
    print(f"Maximum Water Left = {initial_water:.2f} - ({' + '.join(equation_parts)})")
    print(f"Maximum Water Left = {initial_water:.2f} - {total_water_consumed:.2f}")
    print(f"Maximum Water Left = {final_water_left:.2f} liters")
    
    # Required final answer format
    print(f"\n<<<{final_water_left:.2f}>>>")

if __name__ == '__main__':
    # --- User-defined variables ---
    # n: number of 100-liter containers
    # m: distance to destination in km
    n_loads = 3
    m_distance = 80

    solve_horse_problem(n_loads, m_distance)
