import math

def solve_horse_problem(n, m):
    """
    Calculates the maximum amount of water left after a horse travels m km.

    Args:
        n (int): The multiplier for the initial water supply (n * 100 liters).
        m (float): The distance to the destination in kilometers.
    """
    if n * 100 <= m:
        print("Error: The total amount of water is not enough to reach the destination.")
        if n > 0:
            print(f"With {n*100}L of water, the maximum possible distance is less than {n*100} km.")
        return

    # --- Calculation Step ---
    # Find k, the number of depots, and D_k, the distance to the k-th depot.
    k = 0
    D_k = 0.0
    sum_terms = []
    sum_values = []

    # Loop to find which segment the destination 'm' falls into.
    # The loop determines 'k', the number of 100L loads fully utilized for transport.
    while True:
        # The consumption rate for the next segment is (2*(n-k) - 1)
        denominator = 2 * (n - k) - 1
        
        # If denominator is 0 or less, it means we are on the last leg with <= 100L of water.
        if denominator <= 0:
            break
            
        # Calculate the length of the next potential segment
        d_next = 100.0 / denominator
        
        # If m is within this next segment, we've found our k. Break the loop.
        if m < D_k + d_next:
            break
        
        # Otherwise, complete this segment and update k and D_k
        k += 1
        D_k += d_next
        sum_terms.append(f"1/{int(2 * n - 2 * k + 1)}")
        sum_values.append(1.0 / (2 * n - 2 * k + 1))

    # --- Output Step ---
    print("This problem is solved by finding the optimal strategy for transporting water, which involves creating depots.")
    print("The final amount of water is calculated using a formula derived from this strategy.\n")
    
    print("--------------------------------------------------")
    print(f"INPUTS: n = {n} (initial water = {n*100}L), m = {m} km")
    print("--------------------------------------------------\n")

    print("Step 1: Determine 'k', the number of 100L water loads fully consumed to create depots before the final leg.")
    print(f"The calculation finds that k = {k}.\n")

    print(f"Step 2: Calculate D_k, the total distance to the {k}-th depot.")
    print("The formula for D_k is: D_k = 100 * Sum_{i=1 to k} (1 / (2*n - 2*i + 1))")
    
    if k == 0:
        print("Since k=0, no depots were established. D_0 = 0 km.\n")
        D_k_val = 0.0
    else:
        sum_str = " + ".join(sum_terms)
        D_k_val = 100 * sum(sum_values)
        print(f"For k={k}, the summation is: 100 * ( {sum_str} )")
        print(f"D_{k} = {D_k_val:.4f} km.\n")

    print("Step 3: Calculate the maximum water left using the final formula.")
    print("Formula: W_left = 100*(n - k) - (2*(n - k) - 1) * (m - D_k)\n")

    print("Substituting the determined values into the equation:")
    print(f"W_left = 100*({n} - {k}) - (2*({n} - {k}) - 1) * ({m} - {D_k_val:.4f})")
    
    # Breaking down the calculation
    n_minus_k = n - k
    consumption_rate = 2 * n_minus_k - 1
    remaining_dist = m - D_k_val
    
    print(f"W_left = 100*({n_minus_k}) - ({consumption_rate}) * ({remaining_dist:.4f})")
    
    term1 = 100 * n_minus_k
    term2 = consumption_rate * remaining_dist
    water_left = term1 - term2
    
    print(f"W_left = {term1:.4f} - {term2:.4f}")
    
    print("\n--------------------------------------------------")
    print(f"Final Answer: The maximum amount of water left is {water_left:.4f} liters.")
    print("--------------------------------------------------")
    
    global final_answer_value
    final_answer_value = water_left

# --- User Inputs ---
# You can change these values to solve for different scenarios
n_val = 4
m_val = 100

solve_horse_problem(n_val, m_val)

# Return the final answer in the specified format
print(f"\n<<< {final_answer_value:.4f} >>>")