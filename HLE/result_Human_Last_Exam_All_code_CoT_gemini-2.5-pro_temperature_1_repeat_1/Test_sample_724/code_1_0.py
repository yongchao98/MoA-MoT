import math

def solve_horse_problem(n, m):
    """
    Calculates the maximum amount of water left after a horse travels m km.

    Args:
        n: An integer, representing the initial n*100 liters of water.
        m: An integer or float, the distance to the destination in km.
    """
    if n * 100 < m:
        print("The destination is unreachable with the given amount of water.")
        return

    # Find k, the number of trips required for the final leg of the journey.
    # We loop from j=n down to 1. Each j represents a stage where j*100
    # liters are being moved.
    k = 0
    dist_covered_so_far = 0.0
    for j in range(n, 0, -1):
        # Max distance that can be covered in this stage (while moving j*100 liters)
        max_dist_in_stage = 100.0 / (2 * j - 1)
        if m < dist_covered_so_far + max_dist_in_stage:
            k = j
            break
        dist_covered_so_far += max_dist_in_stage
    
    # If loop finishes without break, it means m equals the max possible distance with depots
    if k == 0:
        k = 1 

    # Calculate the summation part of the formula: sum_{j=k+1 to n} (1/(2j-1))
    sum_val = 0.0
    sum_str_parts = []
    for j in range(k + 1, n + 1):
        denominator = 2 * j - 1
        sum_val += 1.0 / denominator
        sum_str_parts.append(f"1/{denominator}")

    sum_str = " + ".join(sum_str_parts) if sum_str_parts else "0"

    # The distance covered before the final leg starts
    dist_before_final_leg = 100.0 * sum_val
    
    # Water at the start of the final leg
    water_at_start_final_leg = k * 100.0
    
    # Distance to travel in the final leg
    dist_in_final_leg = m - dist_before_final_leg
    
    # Consumption rate in the final leg
    consumption_rate = 2 * k - 1
    
    # Water consumed in the final leg
    water_consumed_final_leg = dist_in_final_leg * consumption_rate
    
    # Final amount of water left
    water_left = water_at_start_final_leg - water_consumed_final_leg

    # Construct the equation string for the final output
    equation_str = f"{k}*100 - ({m} - 100 * ({sum_str})) * {consumption_rate}"
    
    print("This problem can be solved with the general formula:")
    print("W_left = k*100 - (m - 100 * \u03A3[j=k+1 to n] 1/(2j-1)) * (2k-1)\n")
    print(f"For n={n} and m={m}, we find k={k}.")
    print("Plugging the numbers into the formula:")
    print(f"Maximum water left = {equation_str}")
    print(f"                     = {water_left:.4f} liters")


# --- User-configurable values ---
# n: initial water is n * 100 liters
n = 4
# m: distance to destination in km
m = 100
# ------------------------------------

solve_horse_problem(n, m)