import math

def solve_horse_problem(n, m):
    """
    Calculates the maximum amount of water left after a horse travels m km.

    Args:
        n (int): The number of 100-liter loads of water available at the origin.
        m (float): The total distance to the destination in km.
    """
    # Assumption check from the problem statement
    # The max distance possible is calculated to ensure the trip is possible.
    max_dist = 0
    if n > 1:
        # Sum of distances for all stages except the last one (which has 1-way trip)
        for i in range(n - 1):
            num_loads = n - i
            consumption_rate = 2 * num_loads - 1
            max_dist += 100.0 / consumption_rate
    # The last 100L can be carried 100 km.
    max_dist += 100.0
    
    if m > max_dist:
        print(f"Error: The destination at {m} km is unreachable with {n*100} liters of water.")
        print(f"The maximum possible distance is approximately {max_dist:.2f} km.")
        return

    print("### Solving the Horse and Water Problem ###\n")
    print("The optimal strategy is to move water forward in stages, creating depots.")
    print("At each stage, we consume exactly 100L of water, reducing the number of round trips required.\n")
    
    # --- Step 1: Find p, the stage in which the destination m is reached ---
    p = 1
    sum_of_distances_before_p = 0.0
    summation_string_parts = []
    
    while True:
        # Number of loads at the beginning of stage p
        num_loads = n - p + 1
        
        # Consumption rate (L/km) for stage p
        consumption_rate = 2 * num_loads - 1
        
        # If consumption rate is 0 or less, it means we have no water to move.
        # This case is handled by the max_dist check above.
        if consumption_rate <= 0:
            break
            
        # Distance covered if we were to complete stage p
        dist_of_stage_p = 100.0 / consumption_rate
        
        if m < sum_of_distances_before_p + dist_of_stage_p:
            # Destination is reached within this stage p. We found p.
            break
        else:
            # This stage is fully completed.
            sum_of_distances_before_p += dist_of_stage_p
            summation_string_parts.append(f"100/({consumption_rate})")
            p += 1

    print(f"Calculation for n={n} (initial water = {n*100}L) and m={m}km:\n")
    print(f"Step 1: Determine the final stage 'p'.")
    print(f"The destination is reached during stage p = {p}.\n")
    
    # --- Step 2: Calculate the water left using the derived formula ---
    
    print("Step 2: Use the formula for remaining water, expressed with a summation.")
    print("Formula: W_left = (n - p + 1)*100 - (m - Sum) * (2*n - 2*p + 1)")
    print("Where 'Sum' is the total distance of the completed stages before stage p:")
    
    if p == 1:
        print("Sum = 0 (since destination is in the first stage)")
    else:
        sum_str = " + ".join(summation_string_parts)
        print(f"Sum = {sum_str}")
        print(f"Sum = {sum_of_distances_before_p:.2f} km")
    
    # Calculate each term in the formula
    final_stage_loads = n - p + 1
    initial_water_final_stage = final_stage_loads * 100
    
    distance_in_final_stage = m - sum_of_distances_before_p
    
    consumption_rate_final_stage = 2 * final_stage_loads - 1
    
    water_consumed_final_stage = distance_in_final_stage * consumption_rate_final_stage
    
    water_left = initial_water_final_stage - water_consumed_final_stage
    
    print("\nStep 3: Plug the numbers into the formula:\n")

    # Output each number in the final equation
    print(f"Water Left = ({n} - {p} + 1)*100 - ({m} - {sum_of_distances_before_p:.2f}) * (2*{n} - 2*{p} + 1)")
    print(f"Water Left = ({final_stage_loads})*100 - ({distance_in_final_stage:.2f}) * ({consumption_rate_final_stage})")
    print(f"Water Left = {initial_water_final_stage} - {water_consumed_final_stage:.2f}")
    print(f"Water Left = {water_left:.2f} liters")
    
    print(f"\nFinal Answer: The maximum amount of water left is {water_left:.2f} liters.")
    
    # Return the final numeric answer for the <<<>>> tag
    return water_left

# Example Usage with values from the prompt's context
# You can change these values to solve for different scenarios
n_loads = 3 # Represents 300 liters of water
m_distance = 85 # Represents 85 km distance

# --- Execute the solver ---
final_result = solve_horse_problem(n=n_loads, m=m_distance)
print(f"\n<<<{final_result:.2f}>>>")
