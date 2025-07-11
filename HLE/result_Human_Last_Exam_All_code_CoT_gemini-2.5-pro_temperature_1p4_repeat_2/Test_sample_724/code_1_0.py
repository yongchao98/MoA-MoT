def solve_horse_water_problem(n, m):
    """
    Calculates the maximum amount of water left after a horse travels m km.

    Args:
        n (int): The number of 100-liter water loads at the start.
        m (float): The total distance to travel in kilometers.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.")
        return
    if not isinstance(m, (int, float)) or m < 0:
        print("Error: m must be a non-negative distance.")
        return
        
    initial_water = float(n * 100)
    
    # Check if the journey is possible in the first place, even with a simple one-way trip.
    # A more rigorous check is done implicitly by the loop.
    if initial_water <= m:
         max_dist = 0
         # Calculate the absolute maximum range with n loads.
         for i in range(n, 0, -1):
             max_dist += 100.0 / (2*i - 1)
         if m > max_dist:
             print(f"Destination is unreachable. The maximum possible distance with {n} loads is {max_dist:.2f} km.")
             return


    water_at_leg_start = initial_water
    distance_traveled = 0.0
    
    # Iterate through the legs of the journey, from k=n loads down to 2 loads
    for k in range(n, 1, -1):
        # Consumption rate is (2k-1) liters per km for this leg
        consumption_rate = 2 * k - 1
        
        # The distance covered by consuming one 100L load at this rate
        leg_distance = 100.0 / consumption_rate
        
        # Check if the destination m falls within this leg
        if m <= distance_traveled + leg_distance:
            distance_in_leg = m - distance_traveled
            water_consumed = distance_in_leg * consumption_rate
            final_water = water_at_leg_start - water_consumed
            
            # --- Construct the equation string for the output ---
            # Part 1: Water at the start of the leg
            water_start_str = f"{k}*100"
            
            # Part 2: Summation of distances of previous legs
            sum_terms = [f"100/(2*{j}-1)" for j in range(n, k, -1)]
            sum_of_distances_str = " + ".join(sum_terms) if sum_terms else "0"
            
            # Part 3: Consumption rate for this leg
            rate_str = f"(2*{k}-1)"
            
            # Assemble the final equation
            equation = f"{water_start_str} - ({m} - ({sum_of_distances_str})) * {rate_str}"
            
            print("The maximum amount of water left is calculated by the following equation:")
            print(f"Final Water = {equation}")
            print(f"Which evaluates to: {final_water:.4f} liters")
            return final_water

        # If m is beyond this leg, update state for the next iteration
        distance_traveled += leg_distance
        water_at_leg_start -= 100.0

    # If the loop completes, it means m is in the final one-way leg (k=1)
    consumption_rate = 1
    distance_in_leg = m - distance_traveled
    
    if water_at_leg_start < distance_in_leg:
        # This case should be caught by the initial check, but as a safeguard:
        print(f"Destination is unreachable. Not enough water for the final leg.")
        return 0.0
        
    water_consumed = distance_in_leg * consumption_rate
    final_water = water_at_leg_start - water_consumed

    # --- Construct the equation string for the final leg ---
    # Part 1: Water at the start is 100L
    water_start_str = "100"
    
    # Part 2: Summation of distances of all previous legs
    sum_terms = [f"100/(2*{j}-1)" for j in range(n, 1, -1)]
    sum_of_distances_str = " + ".join(sum_terms)

    # Assemble the final equation
    equation = f"{water_start_str} - ({m} - ({sum_of_distances_str}))"

    print("The maximum amount of water left is calculated by the following equation:")
    print(f"Final Water = {equation}")
    print(f"Which evaluates to: {final_water:.4f} liters")
    return final_water


# --- Example Usage ---
# You can change these values to match your specific problem
# n: number of 100L water loads.
# m: distance to travel in km.
example_n = 4
example_m = 100

# Execute the function with the example values
final_amount = solve_horse_water_problem(example_n, example_m)

# The final numeric answer for the example above
# This line is for the final answer tag and would not normally be part of the script.
if final_amount is not None:
    print(f"\n<<<Result for n={example_n}, m={example_m}>>>")
    print(f'<<<{final_amount:.4f}>>>')
