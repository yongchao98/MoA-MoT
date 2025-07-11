import math

def solve_horse_water_problem(n, m):
    """
    Calculates the maximum amount of water left after a horse travels a distance m.

    Args:
        n (int): The number of 100-liter loads of water at the origin.
        m (float): The total distance to the destination in kilometers.
    """
    if n * 100 <= m:
        print(f"Warning: The initial water {n*100}L is not greater than the distance {m}km.")
        print("According to the problem assumption (n*100 > m), this case is not considered.")
        print("It's likely impossible to reach the destination.")
        # Under the simplified (but incorrect) assumption, any travel consumes at least 1L/km, so it would fail.
        # Let's show a result of 0, though a negative result is more realistic.
        print("\nMaximum water left:")
        print("0.00")
        return 0.0

    # Find the stage 'k' in which the destination 'm' lies.
    # A stage k begins after k-1 depots have been established.
    distance_to_last_depot = 0.0
    k = 1
    for i in range(1, n + 1):
        # Current number of 100L loads to transport
        num_loads = n - i + 1
        
        # If there's only one load left, we are in the final one-way trip phase
        if num_loads <= 1:
            k = i
            break
            
        consumption_rate = 2 * num_loads - 1
        
        # Distance that can be covered in this stage to consume exactly 100L
        stage_distance = 100.0 / consumption_rate
        
        if distance_to_last_depot + stage_distance >= m:
            k = i
            break
        
        distance_to_last_depot += stage_distance
        k = i + 1

    # --- Build the equation string for the final answer ---
    
    # Part 1: The summation part representing the distance to the (k-1)th depot
    summation_terms = []
    if k > 1:
        for i in range(1, k):
            denominator = 2 * (n - i + 1) - 1
            summation_terms.append(f"1/{denominator}")
        summation_expr = f"100 * ({' + '.join(summation_terms)})"
    else:
        # if k=1, the travel starts from origin (x_0 = 0)
        summation_expr = "0"
        
    # Part 2: The full equation string
    # Formula: W(m) = (n-k+1)*100 - (2*(n-k+1)-1) * (m - x_{k-1})
    final_equation = f"Max water = ({n} - {k} + 1) * 100 - (2 * ({n} - {k} + 1) - 1) * ({m} - {summation_expr})"

    # --- Calculate the final numerical result ---
    water_at_start_of_leg = (n - k + 1) * 100.0
    rate_on_leg = 2 * (n - k + 1) - 1
    distance_of_leg = m - distance_to_last_depot
    
    water_left = water_at_start_of_leg - rate_on_leg * distance_of_leg
    
    # --- Print the final output ---
    print("The amount of water left can be expressed with the following formula, where k is the travel stage containing the destination:")
    print(final_equation)
    
    if water_left < 0:
        print("\nNote: The calculated water left is negative, which means the destination is unreachable with the given amount of water.")

    print("\nThe calculated maximum amount of water left is:")
    print(f"{water_left:.2f}")
    
    return water_left

# Example values for n and m
n_val = 4
m_val = 100

# Solve the problem and print the output
final_answer = solve_horse_water_problem(n_val, m_val)
print(f"\n<<<{final_answer:.2f}>>>")