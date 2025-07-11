def solve_horse_water_problem(n, m):
    """
    Calculates the maximum amount of water left after a horse travels m km.

    Args:
        n: An integer representing the number of 100-liter water batches at the start.
        m: An integer representing the total distance to travel in km.
    """
    # The problem assumes n*100 > m, so travel is always possible.

    # k is the number of full stages completed.
    # A stage is completed when 100L of water is consumed to reduce the number of trips.
    k = 0
    dist_covered_by_stages = 0.0
    
    # This list will store the string representation of each term in the summation
    summation_terms = []
    # This list will store the numerical value of each term
    summation_values = []

    # Loop through the n-1 possible depot-laying stages
    for i in range(1, n):
        # Denominator for the distance calculation in the current stage
        # Number of segments to travel for this stage: 2 * (n - i + 1) - 1
        denominator = 2 * (n - i + 1) - 1
        
        # Optimal distance for this stage
        stage_dist = 100.0 / denominator

        # Check if we can complete this stage before reaching m
        if dist_covered_by_stages + stage_dist <= m:
            dist_covered_by_stages += stage_dist
            summation_terms.append(f"100 / {denominator}")
            summation_values.append(f"{stage_dist:.2f}")
            k = i
        else:
            # m is reached within this stage, so we break the loop
            break

    # After the loop, k holds the number of completed stages.
    # dist_covered_by_stages is the total distance of these k stages, D_k.
    
    print(f"For n={n} (initial water = {n*100}L) and destination distance m={m}km:")
    print("-" * 50)

    # Case 1: m is reached after all n-1 stages are complete (k == n - 1).
    # This means we are in the final one-way trip with 100L of water.
    if k == n - 1:
        # Water left = 100 - (remaining distance)
        # remaining distance = m - dist_covered_by_stages
        water_left = 100.0 - (m - dist_covered_by_stages)
        
        print("The destination is reached in the final one-way leg of the journey.")
        print("The calculation is: Water Left = 100 - m + (Sum of stage distances)")
        
        summation_str = " + ".join(summation_terms)
        print(f"Water Left = 100 - {m} + ({summation_str})")
        
        values_str = ' + '.join(summation_values)
        print(f"Water Left = 100 - {m} + ({values_str})")
        
        print(f"Water Left = 100 - {m} + {dist_covered_by_stages:.2f}")
        print(f"Water Left = {100.0 - m:.2f} + {dist_covered_by_stages:.2f}")

    # Case 2: m is reached during an intermediate stage (stage k+1).
    else:
        # Water at the start of this stage (k+1) is (n-k)*100
        water_at_start_stage = (n - k) * 100.0
        
        # Number of segments for travel in this stage
        num_segments = 2 * (n - k) - 1
        
        # Remaining distance to cover in this stage
        remaining_dist = m - dist_covered_by_stages
        
        # Water consumed in this partial stage
        water_consumed_in_stage = num_segments * remaining_dist
        
        # Final water left
        water_left = water_at_start_stage - water_consumed_in_stage
        
        print(f"The destination is reached during stage {k + 1} (after {k} full stages).")
        print("The calculation is: Water Left = Water at Stage Start - Water Consumed in Stage")
        
        if k > 0:
            summation_str = " + ".join(summation_terms)
            print(f"Water Left = ({n}-{k})*100 - (2*({n}-{k})-1) * ({m} - ({summation_str}))")
            values_str = ' + '.join(summation_values)
            print(f"Water Left = {water_at_start_stage} - {num_segments} * ({m} - ({values_str}))")
            print(f"Water Left = {water_at_start_stage} - {num_segments} * ({m} - {dist_covered_by_stages:.2f})")
        else: # k=0, means m is reached in the very first stage
            print(f"Water Left = {n}*100 - (2*{n}-1) * {m}")
        
        print(f"Water Left = {water_at_start_stage:.2f} - {num_segments} * {remaining_dist:.2f}")
        print(f"Water Left = {water_at_start_stage:.2f} - {water_consumed_in_stage:.2f}")

    print("-" * 50)
    print(f"Maximum water left at destination: {water_left:.2f} liters")
    return water_left

# Example Usage:
# n = number of 100-liter water containers
# m = distance to travel in km
n_val = 4
m_val = 75
final_water = solve_horse_water_problem(n_val, m_val)