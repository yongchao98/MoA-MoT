def solve_horse_problem(n, m):
    """
    Calculates the maximum water left and prints the formula.

    Args:
        n (int): The number of 100-liter containers of water.
        m (int): The distance to the destination in kilometers.
    """
    if not isinstance(n, int) or not isinstance(m, (int, float)) or n <= 0 or m < 0:
        print("Please provide a positive integer for n and a non-negative number for m.")
        return

    # Total initial water should be greater than the distance.
    # The core of the problem is about capacity, but this is a basic sanity check.
    if n * 100 < m:
        print(f"Not enough water to travel {m} km. The horse can go at most {n * 100} km with one load.")
        return

    # Find the final stage k, where we are shuttling k loads.
    # We start with k_iter = n, and find which stage the distance m falls into.
    k_final = 0
    m_remaining = float(m)
    for k_iter in range(n, 0, -1):
        # The journey must end in the final stage with 1 load if m is large enough.
        if k_iter == 1:
            k_final = 1
            break
        
        # Calculate the maximum distance coverable in the current stage (k_iter)
        # by consuming one full 100L load.
        stage_dist = 100.0 / (2 * k_iter - 1)
        
        # If the remaining distance is less than or equal to the distance of this stage,
        # it means the journey ends within this stage.
        if m_remaining <= stage_dist:
            k_final = k_iter
            break
        
        # If not, we complete this stage and move to the next.
        m_remaining -= stage_dist

    # Now, construct the equation string.
    k = k_final
    
    # This is the distance to the start of the final stage k.
    # It's the sum of distances of all prior stages (from n down to k+1).
    depot_dist_parts = []
    # Summation runs from j=n down to k+1.
    for j in range(n, k, -1):
        depot_dist_parts.append(f"100/(2*{j}-1)")
        
    if not depot_dist_parts:
        depot_dist_str = "0"
    else:
        # We wrap the sum in parentheses for clarity if there are multiple terms.
        depot_dist_str = " + ".join(depot_dist_parts)
        if len(depot_dist_parts) > 1:
            depot_dist_str = f"({depot_dist_str})"

    # The final equation is: water_at_start_of_stage - water_consumed_in_stage
    # water_at_start_of_stage = k * 100
    # water_consumed_in_stage = (m - distance_to_depot) * consumption_rate
    final_equation = f"{k}*100 - ({m} - {depot_dist_str}) * (2*{k}-1)"
    
    print("The maximum amount of water left at the destination is given by the following equation:")
    print(final_equation)


# --- Example Usage ---
# You can change the values of n and m here.
# n = number of 100-liter loads
# m = distance in km
example_n = 4
example_m = 100

solve_horse_problem(example_n, example_m)
