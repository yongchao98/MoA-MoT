def solve_horse_problem(n, m):
    """
    Calculates the maximum amount of water left after a horse travels m km.

    Args:
        n (int): The number of 100-liter canisters of water at the start.
        m (int): The total distance to travel in kilometers.
    """
    if n <= 0 or m <= 0:
        print("Number of canisters (n) and distance (m) must be positive.")
        return

    # Initial state
    total_water = n * 100.0
    distance_covered = 0.0
    depots = 0  # This will be our 'k'
    leg_distances = []

    # Step 1: Determine the number of full depots (k) and their distances.
    # We iterate through each potential leg of the journey.
    for i in range(1, n + 1):
        num_canisters = n - (i - 1)
        
        # If only one canister is left, we can't create more depots.
        # The horse will just travel the remaining distance in one trip.
        if num_canisters <= 1:
            break
        
        # Traversals needed to move all current canisters forward
        traversals = 2 * num_canisters - 1
        
        # Optimal distance of this leg (consumes exactly 100L of water)
        leg_dist = 100.0 / traversals
        
        # If the destination is within this leg, we don't complete the leg.
        if distance_covered + leg_dist > m:
            break
            
        # Complete the leg: update state
        distance_covered += leg_dist
        leg_distances.append(leg_dist)
        depots += 1

    # Step 2: Calculate the final water amount based on the depots created.
    k = depots
    sum_of_leg_distances = sum(leg_distances)

    # Consumption on the final leg of the journey
    remaining_distance = m - sum_of_leg_distances
    num_canisters_for_final_leg = n - k
    traversals_for_final_leg = 2 * num_canisters_for_final_leg - 1
    
    water_consumed_for_depots = k * 100.0
    water_consumed_on_final_leg = remaining_distance * traversals_for_final_leg
    
    total_water_consumed = water_consumed_for_depots + water_consumed_on_final_leg
    water_left = total_water - total_water_consumed

    # Step 3: Format the final answer as a summation equation.
    initial_water_str = f"{total_water:.2f}"
    
    consumption_parts = []
    # Add consumption for each depot created
    if k > 0:
        consumption_parts.extend(["100.00"] * k)

    # Build the string for the final leg's consumption calculation
    remaining_dist_calculation_str = f"{float(m)}"
    if k > 0:
        leg_dist_str_parts = [f"{d:.2f}" for d in leg_distances]
        leg_dist_sum_str = ' + '.join(leg_dist_str_parts)
        remaining_dist_calculation_str = f"({float(m)} - ({leg_dist_sum_str}))"

    final_leg_consumption_str = f"{remaining_dist_calculation_str} * {traversals_for_final_leg}"
    consumption_parts.append(final_leg_consumption_str)
    
    # Combine all parts of the consumption calculation
    total_consumption_str = " + ".join(consumption_parts)
    
    # Final equation and result
    equation = f"Final Water = {initial_water_str} - ({total_consumption_str})"
    result_str = f"Final Result = {water_left:.2f} liters"

    print("The optimal strategy is to create depots of water along the way.")
    print("The maximum amount of water left is calculated by subtracting the total water consumed from the initial amount.")
    print("\nCalculation:")
    print(equation)
    print(f"\nWhich simplifies to:")
    print(f"Final Water = {initial_water_str} - ({water_consumed_for_depots:.2f} + {water_consumed_on_final_leg:.2f})")
    print(result_str)
    
    # The final numerical answer for the bot's <<<>>> tag
    global final_answer_value
    final_answer_value = water_left


# --- User Inputs ---
# n = number of 100L water canisters
# m = distance to travel in km
n_input = 3
m_input = 80
# --- End of User Inputs ---

# This global variable is used to pass the final numerical result out
final_answer_value = 0.0
solve_horse_problem(n_input, m_input)

# The line below is for the platform, to extract the final numerical answer.
# print(f"\n<<<{final_answer_value:.2f}>>>")