import math

def solve_horse_problem(n, m):
    """
    Calculates the maximum amount of water left after a horse travels m km.
    Also prints the step-by-step formula as requested.

    Args:
        n (int): The number of 100-liter water units available at the origin.
        m (float): The distance to travel in kilometers.
    """

    # First, calculate the absolute maximum distance possible to check if m is reachable.
    max_dist = 0
    phase_dists = []
    for i in range(1, n + 1):
        # In phase 'i', k = n - i + 1 units of water are being moved.
        k = n - i + 1
        consumption_rate = 2 * k - 1
        # The distance covered in this phase consumes exactly 100L of water.
        dist = 100.0 / consumption_rate if consumption_rate > 0 else float('inf')
        phase_dists.append(dist)
        max_dist += dist

    if m > max_dist:
        print(f"Error: The destination m={m} km is unreachable with {n*100}L of water.")
        print(f"The maximum possible distance is {max_dist:.2f} km.")
        return 0.0

    # Step 1: Find the phase 'p' corresponding to the distance 'm'.
    cumulative_dist = 0.0
    final_p = 0
    for i in range(n):
        p_loop = i + 1
        current_phase_dist = phase_dists[i]
        if m < cumulative_dist + current_phase_dist:
            final_p = p_loop
            break
        cumulative_dist += current_phase_dist

    # Handle the case where m lands exactly at the end of the last phase
    if final_p == 0 and math.isclose(m, max_dist):
        final_p = n

    # Step 2: Set up the variables for the formula based on the phase 'final_p'.
    k = n - final_p + 1  # Number of 100L units being transported in this phase
    consumption_rate = 2 * k - 1
    
    # Step 3: Construct the summation part of the equation for printing.
    # The summation represents the total distance to the start of the current phase.
    sum_val = 0.0
    sum_str_parts = []
    sum_calc_parts = []

    for i in range(1, final_p):
        k_sum = n - i + 1
        C_sum = 2 * k_sum - 1
        sum_str_parts.append(f"100/(2*{k_sum}-1)")
        sum_calc_parts.append(f"100/{C_sum}")
        sum_val += 100.0 / C_sum

    if not sum_str_parts:
        sum_str_expression = "0"
    else:
        sum_str_expression = f"({' + '.join(sum_str_parts)})"

    # Step 4: Calculate the final water amount.
    water_left = (k * 100.0) - (m - sum_val) * consumption_rate

    # Step 5: Print the explanation and the final equation with all numbers.
    print(f"For n = {n} (initial water = {n*100}L) and distance m = {m}km:")
    print("The amount of water left (W) is given by the formula:")
    print("W = (Water at start of phase) - (Water consumed in phase)")
    
    print("\n--- Equation in summation form ---")
    if final_p > 1:
        summation_symbolic = f"sum_{{i=1}}^{{{final_p - 1}}} (100 / (2*({n}-i+1)-1))"
        print(f"W = ({k} * 100) - (m - {summation_symbolic}) * (2*{k}-1)")
    else: # In phase 1, the summation term is zero
        print(f"W = ({k} * 100) - m * (2*{k}-1)")

    print("\n--- Equation with numbers plugged in ---")
    print(f"W = ({k} * 100) - ({m} - {sum_str_expression}) * {consumption_rate}")

    print("\n--- Calculation Steps ---")
    dist_at_phase_start_str = "0"
    if sum_calc_parts:
         dist_at_phase_start_str = f"{' + '.join(sum_calc_parts)}"
    
    print(f"1. Distance to start of phase {final_p}: {dist_at_phase_start_str} = {sum_val:.2f} km")
    water_at_phase_start = k * 100
    print(f"2. Water at start of phase {final_p}: {water_at_phase_start} L")
    dist_in_phase = m - sum_val
    print(f"3. Distance traveled within phase {final_p}: {m} - {sum_val:.2f} = {dist_in_phase:.2f} km")
    consumption_in_phase = dist_in_phase * consumption_rate
    print(f"4. Water consumed in phase: {dist_in_phase:.2f} km * {consumption_rate} L/km = {consumption_in_phase:.2f} L")
    print("\n--- Final Calculation ---")
    print(f"W = {water_at_phase_start} - {consumption_in_phase:.2f}")
    print(f"W = {water_left:.2f} liters")

    return water_left

if __name__ == '__main__':
    # --- You can change these values ---
    # n: number of 100-liter water portions. e.g., n=3 means 300 liters.
    # m: distance to the destination in kilometers.
    n = 4
    m = 100
    
    final_water = solve_horse_problem(n, m)
