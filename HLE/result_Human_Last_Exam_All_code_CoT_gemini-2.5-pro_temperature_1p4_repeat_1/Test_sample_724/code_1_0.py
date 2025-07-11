import math

def solve_horse_water_problem(n, m):
    """
    Calculates the maximum amount of water left after a horse trip.

    Args:
        n: The number of 100-liter loads of water at the origin.
        m: The total distance to the destination in kilometers.
    """
    if n * 100 < m:
        print("The horse does not have enough water to complete the journey.")
        return

    # Determine the final stage 'j'
    m_rem = float(m)
    full_stages_dist_sum = 0.0
    # A list to hold the terms for the summation string
    full_stages_sum_terms = []
    
    j = 1 # Assume final leg is with 1 load by default

    for k in range(n, 1, -1):
        # consumption_rate is 2k-1 for k > 1 loads
        consumption_rate = 2 * k - 1
        # The optimal distance for the stage where we move k loads
        stage_dist = 100.0 / consumption_rate
        
        if m_rem > stage_dist:
            m_rem -= stage_dist
            full_stages_dist_sum += stage_dist
            # Store the term for the equation string, like "100/(2*3-1)"
            full_stages_sum_terms.append(f"100/(2*{k}-1)")
        else:
            # The journey ends in the stage with k loads
            j = k
            break
            
    print(f"The horse starts with {n*100} liters of water to travel {m} km.")
    print(f"The journey will end during the stage where the horse is moving {j} x 100L loads of water.\n")

    # Now, build the equation and calculate the result based on 'j'
    
    if j > 1:
        # Case 1: Final stage has j > 1 loads
        consumption_rate = 2 * j - 1
        
        # Build the summation part of the equation string
        if not full_stages_sum_terms:
            sum_str = "0"
        else:
            sum_str = " + ".join(full_stages_sum_terms)
        
        # Construct the final equation string
        equation = f"Water Left = {j}*100 - (2*{j}-1) * ({m} - ({sum_str}))"
        
        # Calculate the result
        water_left = j * 100.0 - consumption_rate * (m - full_stages_dist_sum)
        
    else: # j == 1
        # Case 2: Final stage has j = 1 load
        sum_terms = []
        total_sum_val = 0.0
        for i in range(1, n + 1):
            term_val = 1.0 / (2 * i - 1)
            total_sum_val += term_val
            sum_terms.append(f"1/(2*{i}-1)")

        sum_str = " + ".join(sum_terms)
        
        # Construct the final equation string
        equation = f"Water Left = 100 * ({sum_str}) - {m}"
        
        # Calculate the result
        water_left = 100.0 * total_sum_val - m

    print("The amount of water left is calculated by the following equation:")
    print(equation)
    print(f"\nMaximum amount of water left at the destination: {water_left:.2f} liters")
    
# --- User inputs ---
# n: number of 100L water loads.
# m: distance in kilometers.
n = 4
m = 200

# Solve the problem with the given inputs
solve_horse_water_problem(n, m)
<<<25.93>>>