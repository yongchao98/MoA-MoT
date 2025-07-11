import math

def calculate_max_water(n, m):
    """
    Calculates the maximum amount of water left after a horse travels m km.

    Args:
        n (int): The number of 100-liter containers of water initially available.
        m (float): The total distance to travel in kilometers.
    """
    # Assumption check from the problem statement
    if n * 100 <= m:
        print("It's not possible to reach the destination as the total water is less than or equal to the distance.")
        # Maximum possible travel distance to check for reachability
        max_dist = 0
        for i in range(1, n):
            max_dist += 100.0 / (2 * (n - (i - 1)) - 1)
        max_dist += 100.0 # From the last container
        if m > max_dist:
             print(f"The destination is unreachable. Maximum possible distance is {max_dist:.2f} km.")
             return

    # Step 1: Find k, the number of full segments where a 100-liter cache is dropped.
    k = 0
    # D_k is the total distance covered by these k segments
    D_k = 0.0
    for i in range(1, n):
        # Denominator is 2 * (number of containers being moved) - 1
        # In segment i, we are moving n-(i-1) containers
        segment_dist = 100.0 / (2 * n - 2 * i + 1)
        if D_k + segment_dist < m:
            D_k += segment_dist
            k += 1
        else:
            break

    # Step 2: Calculate the final result using the formula
    rem_containers = n - k
    consumption_rate_final_leg = 2 * rem_containers - 1 if rem_containers > 0 else 1.0
    final_water_amount = 100.0 * rem_containers - (m - D_k) * consumption_rate_final_leg

    # Step 3: Construct the detailed equation string as requested
    equation_str = ""
    if k == 0:
        # This is the case where the destination is within the first segment
        equation_str = f"100 * {n} - {m} * (2 * {n} - 1) = {final_water_amount:.4f}"
    else:
        # Build the string for the summation part D_k = Sum_{i=1 to k} ...
        sum_terms_list = [f"100/(2*{n} - 2*{i} + 1)" for i in range(1, k + 1)]
        sum_str = " + ".join(sum_terms_list)

        term1_str = f"100*({n} - {k})"
        # Add extra parentheses for sums with more than one term for clarity
        term2_str = f"({m} - ({sum_str}))" if k > 1 else f"({m} - {sum_str})"
        term3_str = f"(2*({n} - {k}) - 1)"
        
        equation_str = f"{term1_str} - {term2_str} * {term3_str} = {final_water_amount:.4f}"
        
    print("The final amount of water is calculated by the following equation:")
    print(equation_str)


# --- User-defined values ---
# You can change these values to fit your problem.
# n: represents the initial amount of water in units of 100 liters (e.g., n=3 means 300 liters)
# m: the distance to the destination in kilometers.
n_val = 3
m_val = 50

# Execute the calculation and print the result
calculate_max_water(n_val, m_val)
