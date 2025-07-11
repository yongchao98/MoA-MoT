import math

def solve_horse_problem(n, m):
    """
    Calculates the maximum amount of water left after a horse travels m km.

    Args:
        n (int): The number of 100-liter water caches available at the start.
        m (float): The total distance to the destination in kilometers.
    """
    if n * 100 < m:
        print("Not enough water to reach the destination.")
        return

    # These variables track the state as we complete full stages
    dist_covered_full_stages = 0.0
    stages_completed = 0
    # This list will store the terms for the summation string
    dist_terms_str = []
    
    # This will be the consumption rate for the final leg of the journey
    final_consumption_rate = 0

    # This flag helps determine if the destination is reached
    # while multiple trips are still required.
    destination_reached_in_multitrip = False

    # Iterate backwards from n loads down to 2 loads.
    # The last leg (from 100L down to 0) is handled separately.
    for k in range(n, 1, -1):
        # Consumption rate for moving k caches of 100L water
        consumption_rate = 2 * k - 1
        
        # The distance we can travel in this stage before consuming 100L of water
        leg_dist = 100.0 / consumption_rate
        
        # Check if the destination is within this leg
        if dist_covered_full_stages + leg_dist >= m:
            final_consumption_rate = consumption_rate
            destination_reached_in_multitrip = True
            break
        else:
            # Complete the full leg
            dist_covered_full_stages += leg_dist
            stages_completed += 1
            dist_terms_str.append(f"100/{consumption_rate}")

    # If the loop finished, it means the destination is in the final leg (k=1)
    if not destination_reached_in_multitrip:
        final_consumption_rate = 1.0

    # Calculate the water available at the beginning of the final leg
    water_at_start_of_final_leg = (n - stages_completed) * 100.0
    
    # Calculate the remaining distance to travel
    remaining_distance = m - dist_covered_full_stages
    
    # Calculate the water consumed in this final leg
    water_consumed_in_final_leg = remaining_distance * final_consumption_rate
    
    # The final amount of water left
    final_water_left = water_at_start_of_final_leg - water_consumed_in_final_leg
    
    # Prepare the summation string for the equation
    dist_sum_str = " + ".join(dist_terms_str) if dist_terms_str else "0"

    # Print the final equation as requested
    print(f"For n={n} and m={m}:")
    # Using 'g' format specifier to avoid unnecessary trailing zeros.
    print(f"Maximum amount of water left = {water_at_start_of_final_leg:g} - ({m:g} - ({dist_sum_str})) * {final_consumption_rate:g} = {final_water_left:.4f}")
    
    return final_water_left

# --- User-configurable inputs ---
# n: an integer, representing n*100 liters of water
n = 3
# m: a number, representing the distance to the destination in km
m = 100
# ----------------------------------

final_answer = solve_horse_problem(n, m)
print(f"<<<{final_answer:.4f}>>>")