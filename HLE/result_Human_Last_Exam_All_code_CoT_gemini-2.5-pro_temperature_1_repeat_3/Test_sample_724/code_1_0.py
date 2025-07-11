def solve_horse_problem(n, m):
    """
    Calculates the maximum amount of water left after a journey.

    Args:
        n (int): The number of 100-liter canisters of water at the origin.
        m (float): The total distance to the destination in kilometers.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.")
        return
    if not isinstance(m, (int, float)) or m < 0:
        print("Error: m must be a non-negative number.")
        return
    if n * 100 < m:
        print("Error: Not enough water to complete the journey even in a single trip.")
        # The problem statement assumes n*100 > m, but this is a good check.
        return

    # --- Value Calculation Variables ---
    initial_water_val = n * 100.0
    dist_covered_val = 0.0
    consumed_stages_val = [] # Stores the numeric value of water consumed in each full stage

    # --- String Representation Variables ---
    consumed_stages_str = [] # Stores the string "100" for each full stage
    dist_stages_str = []     # Stores the string representation of each full stage's distance

    # Loop through the multi-canister stages, from n canisters down to 2
    for k in range(n, 1, -1):
        num_trips = 2 * k - 1
        # Max distance for this stage to consume exactly 100L of water
        stage_max_dist_val = 100.0 / num_trips
        
        remaining_dist_to_destination = m - dist_covered_val

        if remaining_dist_to_destination > stage_max_dist_val:
            # We travel the full length of this optimal stage
            dist_covered_val += stage_max_dist_val
            consumed_stages_val.append(100.0)
            consumed_stages_str.append("100")
            dist_stages_str.append(f"100/{num_trips}")
        else:
            # The journey ends in this multi-trip stage.
            # Calculate final water amount (value)
            consumed_in_this_stage_val = num_trips * remaining_dist_to_destination
            total_consumed_val = sum(consumed_stages_val) + consumed_in_this_stage_val
            final_water_val = initial_water_val - total_consumed_val

            # Construct the equation string
            sum_prev_dist_str = " + ".join(dist_stages_str) if dist_stages_str else "0"
            consumed_in_this_stage_str = f"{num_trips} * ({m} - ({sum_prev_dist_str}))"
            
            all_consumed_parts_str = consumed_stages_str + [consumed_in_this_stage_str]
            total_consumption_str = " + ".join(all_consumed_parts_str)
            
            equation_str = f"Maximum water left = {n}*100 - ({total_consumption_str})"
            
            print(equation_str)
            print(f"= {final_water_val:.4f}")
            return

    # If the loop completes, the journey ends in the final 1-canister, one-way trip stage.
    # Calculate final water amount (value)
    final_leg_dist_val = m - dist_covered_val
    total_consumed_val = sum(consumed_stages_val) + final_leg_dist_val
    final_water_val = initial_water_val - total_consumed_val

    # Construct the equation string
    sum_prev_dist_str = " + ".join(dist_stages_str) if dist_stages_str else "0"
    consumed_in_final_leg_str = f"({m} - ({sum_prev_dist_str}))"
    
    all_consumed_parts_str = consumed_stages_str + [consumed_in_final_leg_str]
    total_consumption_str = " + ".join(all_consumed_parts_str)

    equation_str = f"Maximum water left = {n}*100 - ({total_consumption_str})"

    print(equation_str)
    print(f"= {final_water_val:.4f}")


# --- USER INPUT ---
# Set the number of 100-liter water canisters (n) and distance (m) here.
n = 4  # e.g., 4 * 100 = 400 liters of water
m = 100  # e.g., 100 km distance

solve_horse_problem(n, m)
>>> Maximum water left = 4*100 - (100 + 100 + 100 + (100 - (100/7 + 100/5 + 100/3)))
= 70.9524