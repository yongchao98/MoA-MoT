def solve_horse_water_problem(n, m):
    """
    Calculates the maximum amount of water left after a horse travels m km.

    Args:
        n: The number of 100-liter units of water available at the origin.
        m: The total distance to the destination in kilometers.
    """
    if n * 100 <= m:
        print("The horse does not have enough water to reach the destination.")
        return 0

    # Step 1: Determine the number of depots and their locations (distances).
    depots_built = 0
    depot_distance_values = []
    m_remaining_for_depots = float(m)

    for i in range(n - 1):  # Can build at most n-1 depots
        current_loads = n - i
        trips = 2 * current_loads - 1
        
        # The optimal distance to place the next depot
        segment_distance = 100.0 / trips
        
        if m_remaining_for_depots > segment_distance:
            # We build a full depot
            depots_built += 1
            depot_distance_values.append(segment_distance)
            m_remaining_for_depots -= segment_distance
        else:
            # The remaining journey is the final leg
            break

    # Step 2: Calculate the total water consumed and the final amount left.
    # Water consumed to create the depots (100L each)
    water_consumed_for_depots = depots_built * 100.0

    # Water consumed on the final leg of the journey
    loads_for_final_leg = n - depots_built
    trips_for_final_leg = 2 * loads_for_final_leg - 1
    distance_of_final_leg = float(m) - sum(depot_distance_values)
    water_consumed_on_final_leg = distance_of_final_leg * trips_for_final_leg

    total_water_consumed = water_consumed_for_depots + water_consumed_on_final_leg
    water_left = (n * 100.0) - total_water_consumed

    # Step 3: Construct the equation string for the output.
    # Format: Water Left = Initial Water - Water for Depots - Water for Final Leg
    
    # Start with initial water
    equation_str = f"{water_left:.4f} = {n}*100"
    
    # Subtract water for depots
    if depots_built > 0:
        equation_str += f" - {depots_built}*100"
        
    # Subtract water for the final leg, expressed in terms of summations
    final_leg_sum_parts = []
    for i in range(depots_built):
        loads = n - i
        final_leg_sum_parts.append(f"100/(2*{loads}-1)")
    
    if not final_leg_sum_parts:
        final_dist_expr = f"{m}"
    else:
        # We show the sum of distances to the last depot
        final_dist_expr = f"({m} - ({' + '.join(final_leg_sum_parts)}))"

    num_final_loads = n - depots_built
    final_trips_expr = f"(2*{num_final_loads}-1)"

    equation_str += f" - {final_dist_expr} * {final_trips_expr}"

    print("The final amount of water is calculated by the following summation:")
    print(equation_str)
    
    return water_left

# --- Example Usage ---
# You can change these values to test different scenarios
n_loads = 3
m_distance = 120

print(f"Solving for n={n_loads} loads and m={m_distance} km distance:")
final_answer = solve_horse_water_problem(n_loads, m_distance)

# Another example for a shorter distance
# print("\n" + "-"*20)
# n_loads = 3
# m_distance = 15
# print(f"Solving for n={n_loads} loads and m={m_distance} km distance:")
# solve_horse_water_problem(n_loads, m_distance)

<<<33.3333>>>