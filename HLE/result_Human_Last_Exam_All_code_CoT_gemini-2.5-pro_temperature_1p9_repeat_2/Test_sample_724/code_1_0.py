import math

def solve_horse_water_problem(n, m):
    """
    Calculates the maximum amount of water left after a horse travels m km.
    It prints the final calculation steps and the result.

    Args:
        n (int): The initial number of 100-liter water loads.
        m (float): The total distance to travel in kilometers.
    """
    
    # Step 1: Find k, the number of loads being moved in the segment where m falls,
    # and calculate x_k, the cumulative distance to the start of that segment.

    k = 0
    # x_k_dist corresponds to sum_{j=k+1 to n} 100 / (2j - 1)
    x_k_dist = 0.0  
    
    cumulative_dist = 0.0
    found_k = False
    
    # Iterate j from n down to 2. 'j' represents the number of 100L loads
    # being transported in a given stage.
    for j in range(n, 1, -1):
        # The length of the stage where 'j' loads are moved until 100L are consumed
        segment_length = 100.0 / (2 * j - 1)
        
        # Check if the destination 'm' falls within this stage
        if m <= cumulative_dist + segment_length:
            k = j
            x_k_dist = cumulative_dist
            found_k = True
            break
            
        cumulative_dist += segment_length
        
    if not found_k:
        # If the loop completes, 'm' is beyond the point where the water is
        # reduced to 100L. The final leg is made with a single load.
        k = 1
        # The starting point for this final leg is at distance x_1,
        # which is the total distance covered by all depot-laying stages.
        x_k_dist = cumulative_dist

    # Step 2: Calculate the final amount of water left using the derived formula:
    # Water Left = (Water at start of final stage) - (Water consumed in final stage)
    # Water Left = k*100 - (m - x_k) * (2k - 1)
    
    water_left = k * 100.0 - (m - x_k_dist) * (2 * k - 1)

    # Step 3: Construct the equation string to display the calculation as a summation.
    
    # Build the string for the summation part: sum_{j=k+1 to n} 100 / (2j - 1)
    summation_terms = []
    # The loop for the summation is from j = n down to k+1
    for j in range(n, k, -1):
        denominator = 2 * j - 1
        term_str = f"100/{denominator}"
        summation_terms.append(term_str)
        
    summation_str = "0"
    if summation_terms:
        # Join the terms with '+' and wrap in parentheses
        summation_str = f"({' + '.join(summation_terms)})"

    # Assemble the final equation string with all numbers evaluated.
    k_term = f"{k}*100"
    m_term = str(m)
    rate_term = str(2*k - 1)

    # The format is: k*100 - (m - SUM) * (2k-1) = RESULT
    print(f"For n={n} loads and m={m} km distance:")
    print("The final equation representing the maximum water left is:")
    print(f"{k_term} - ({m_term} - {summation_str}) * {rate_term} = {water_left:.4f}")
    print(f"\nThe maximum amount of water left is {water_left:.4f} liters.")


# --- User-defined values ---
# n represents the number of 100-liter containers, so initial water is n * 100
n = 4
# m is the distance to the destination in kilometers
m = 150.0
# ---------------------------

solve_horse_water_problem(n, m)
<<<17.6190>>>