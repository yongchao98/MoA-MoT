import math

def calculate_water_left(n, m):
    """
    Calculates the maximum amount of water left after a horse travels m km.

    Args:
        n: The number of 100-liter water containers available initially.
        m: The total distance to travel in kilometers.
    """

    # Check if the journey is possible.
    # The absolute maximum distance is when all water is optimally used.
    # max_dist = 100 * (1/1 + 1/3 + 1/5 + ... + 1/(2n-1))
    max_dist = 100.0 * sum(1.0 / (2 * i - 1) for i in range(1, n + 1))
    if m > max_dist:
        print(f"Error: The destination m={m} km is unreachable.")
        print(f"With n={n}, the maximum possible distance is {max_dist:.2f} km.")
        return

    # Find k, the number of loads at the start of the final travel segment, and D_k.
    k = 0
    distance_before_stage_k = 0.0
    summation_terms_str = []
    
    # We start with n loads. We determine which stage (k) the distance m falls into.
    for current_loads in range(n, 0, -1):
        # Distance covered to consume 100L and reduce loads from k to k-1
        dist_this_stage = 100.0 / (2 * current_loads - 1)
        
        # If m is beyond this stage, we complete it.
        if m > dist_this_stage:
            m -= dist_this_stage
            distance_before_stage_k += dist_this_stage
            summation_terms_str.append(f"100/(2*{current_loads}-1)")
        else:
            # m is within this stage. This is our final leg.
            k = current_loads
            break
            
    # At this point, `k` is the number of loads for the final leg,
    # `distance_before_stage_k` is D_k, and `m` is the remaining distance (m - D_k).
    
    # Calculate final water left using the formula: W = k*100 - (2k-1)*(m_rem)
    water_left = k * 100.0 - (2 * k - 1) * m
    
    # --- Output Generation ---
    print(f"For n={n} and m={m}, the optimal strategy leads to a final stage with k={k} loads.")
    print("The formula for water left is: Water = k*100 - (2*k-1) * (m - D_k)")
    print(f"where D_k is the distance covered before the final stage, D_k = sum_{{j={k+1} to {n}}}(100/(2j-1))")
    
    # Build the string representation for the D_k summation
    if not summation_terms_str:
        Dk_sum_str = "0"
    else:
        # The terms were added for j=n, n-1, ... k+1. Reverse for correct order.
        Dk_sum_str = " + ".join(reversed(summation_terms_str))

    # Output the final equation with all numbers plugged in.
    # We use the original 'm' value for the equation printout.
    original_m = distance_before_stage_k + m

    print("\nFinal equation with substituted values:")
    final_eq_str = (f"Water Left = {k}*100 - (2*{k}-1) * ({original_m} - ({Dk_sum_str}))")
    print(final_eq_str)
    
    print("\nSimplified calculation:")
    calc_str = (f"Water Left = {k*100} - {2*k-1} * ({original_m} - {distance_before_stage_k:.2f}) "
                f"= {k*100} - {2*k-1} * {m:.2f} = {water_left:.2f}")
    print(calc_str)
    
    print(f"\nFinal Answer: The maximum amount of water left is {water_left:.2f} liters.")
    
    # Required final answer format
    print(f"\n<<<{water_left:.2f}>>>")

# --- Example Usage ---
# You can change these values to test different scenarios
n_val = 3   # n*100 liters of water
m_val = 80  # distance to travel in km

calculate_water_left(n_val, m_val)