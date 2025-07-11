import sys

def solve_water_problem(n, m):
    """
    Calculates the maximum amount of water left after a horse travels m km.
    The solution is presented in terms of a summation of water consumed in stages.

    Args:
        n (int): The number of 100-liter units of water available at the start.
        m (float): The total distance to the destination in kilometers.
    """
    # The problem statement assumes n*100 > m, so there's enough water for a simple one-way trip.
    initial_water = float(n * 100)
    distance_to_travel = float(m)

    if initial_water < distance_to_travel:
        print(f"Error: Initial water ({initial_water:.2f}L) is not enough to travel {distance_to_travel:.2f}km even in a single trip.")
        return 0.0

    distance_traveled = 0.0
    tanks_to_move = n
    consumed_parts_str = []
    total_consumed = 0.0

    while distance_traveled < distance_to_travel:
        # Prevent infinite loops on edge cases or if something goes wrong.
        if tanks_to_move <= 0:
            print("Error: Ran out of tanks unexpectedly.")
            break

        # The consumption rate (L/km) is (2 * k - 1), where k is the number of tanks to move.
        consumption_rate = 2 * tanks_to_move - 1
        
        remaining_dist_to_destination = distance_to_travel - distance_traveled
        
        dist_this_leg = 0.0

        if tanks_to_move > 1:
            # For multiple tanks, the most efficient stage distance consumes exactly 100L.
            max_dist_for_full_stage = 100.0 / consumption_rate

            if max_dist_for_full_stage < remaining_dist_to_destination:
                # We travel a full stage, consuming exactly 100L of water.
                dist_this_leg = max_dist_for_full_stage
                tanks_to_move -= 1 # One tank has been used for transport.
            else:
                # The remaining distance is shorter than a full stage. This is the final leg.
                dist_this_leg = remaining_dist_to_destination
        else: # tanks_to_move == 1
            # Only one tank left. This is a simple one-way trip from the current location.
            # Consumption rate is 1 L/km.
            dist_this_leg = remaining_dist_to_destination

        # Stop if no progress is made (avoids issues with floating point precision).
        if dist_this_leg <= 1e-9:
            break

        consumed_this_leg = dist_this_leg * consumption_rate
        total_consumed += consumed_this_leg
        distance_traveled += dist_this_leg
        
        # Build the part of the equation that shows how consumption was calculated for this stage.
        # This fulfills the requirement to "output each number in the final equation".
        consumed_parts_str.append(f"({dist_this_leg:.2f} * {consumption_rate})")
    
    water_left = initial_water - total_consumed

    # Create the final summation output string.
    summation_str = " + ".join(consumed_parts_str)
    
    print("\n--- Calculation Breakdown ---")
    print(f"Initial Water: {initial_water:.2f} L")
    print(f"Destination: {distance_to_travel:.2f} km")
    print("\nWater consumption is calculated in stages. Each term in the sum represents one stage:")
    print("Consumed Water = (distance of stage 1 * rate 1) + (distance of stage 2 * rate 2) + ...")
    
    # Print the final equation in terms of the summation
    print(f"\nFinal Equation:\n{initial_water:.2f} - ({summation_str}) = {water_left:.2f}")

    print(f"\nMaximum amount of water left: {water_left:.2f} liters")
    
    return water_left

if __name__ == '__main__':
    # Set default values for n and m
    n_tanks = 3     # Corresponds to 300 liters
    m_distance = 100.0 # Corresponds to 100 km

    # Allow user to override defaults via command-line arguments
    try:
        if len(sys.argv) == 3:
            n_tanks = int(sys.argv[1])
            m_distance = float(sys.argv[2])
        
        # Execute the calculation
        final_water_amount = solve_water_problem(n_tanks, m_distance)
        print(f"\n<<<%.2f>>>" % final_water_amount)

    except (IndexError, ValueError):
        print("\nUsage: python your_script.py <n_tanks> <m_distance>")
        print(f"Example: python {sys.argv[0]} 3 100")
        print(f"(Running with default values: n={n_tanks}, m={m_distance})")
        final_water_amount = solve_water_problem(n_tanks, m_distance)
        print(f"\n<<<%.2f>>>" % final_water_amount)