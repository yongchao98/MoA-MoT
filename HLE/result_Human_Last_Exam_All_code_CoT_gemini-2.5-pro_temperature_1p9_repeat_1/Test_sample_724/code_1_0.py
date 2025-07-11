import math

def solve_horse_problem(n, m):
    """
    Calculates the maximum amount of water left after a horse travels m km.

    Args:
        n (int): The number of 100-liter buckets of water at the start.
        m (float): The total distance to travel in kilometers.
    """

    print(f"Calculating for n = {n} (initial water = {n*100}L) and distance m = {m}km.")
    print("-" * 40)
    
    # Check if the journey is possible.
    # The maximum possible distance is when all water is consumed.
    max_dist = 0
    for i in range(1, n + 1):
        max_dist += 100 / (2 * i - 1)
    
    if m > max_dist:
        print("The destination is too far. The horse cannot reach it with the given water.")
        print(f"Maximum reachable distance is {max_dist:.2f} km.")
        return

    # --- Step 1: Find j ---
    # Find 'j', the number of trips needed for the final leg of the journey.
    # We do this by summing the max distance coverable in each stage,
    # starting from stage 'n', until the total distance exceeds 'm'.
    j = 0
    dist_covered_before_stage = 0
    for k in range(n, 0, -1):
        # The max distance covered while consuming 100L of water in stage 'k'
        dist_in_stage_k = 100 / (2 * k - 1)
        if m <= dist_covered_before_stage + dist_in_stage_k:
            j = k
            break
        dist_covered_before_stage += dist_in_stage_k
    
    # Handle the case where m is 0
    if m == 0:
      j = n

    # --- Step 2: Build the equation string ---
    # The formula is: Water_Left = j*100 - (2*j-1) * (m - sum_{i=j+1 to n} [100/(2i-1)])
    # We build the summation part of the string first.
    sum_dist_before_j = 0
    sum_str_parts = []
    # The sum goes from i = j+1 up to n.
    for i in range(j + 1, n + 1):
        term_value = 100 / (2 * i - 1)
        sum_dist_before_j += term_value
        sum_str_parts.append(f"100/(2*{i}-1)")

    if not sum_str_parts:
        sum_str = "0"
    else:
        # We reverse because we calculated from k=n downwards but sum is written j+1 to n
        sum_str = " + ".join(sum_str_parts)

    # --- Step 3: Calculate the result and print ---
    water_left = j * 100 - (2 * j - 1) * (m - sum_dist_before_j)

    # Assemble the final, human-readable equation with all numbers substituted.
    final_equation_str = f"{j} * 100 - (2*{j}-1) * ({m} - ({sum_str}))"

    print("The amount of water remaining is expressed by the following equation:")
    print(final_equation_str)
    
    print("\nWhich calculates to:")
    final_result_str = f"Result = {j*100} - ({2*j-1}) * ({m} - {sum_dist_before_j:.2f}) = {water_left:.2f} liters"
    print(final_result_str)


# --- Example Usage ---
# You can change these values to solve for different scenarios.

# n: Represents the initial amount of water in units of 100 liters.
# (e.g., n=3 means 300 liters)
n_val = 4

# m: The distance to the destination in kilometers.
m_val = 80

solve_horse_problem(n=n_val, m=m_val)