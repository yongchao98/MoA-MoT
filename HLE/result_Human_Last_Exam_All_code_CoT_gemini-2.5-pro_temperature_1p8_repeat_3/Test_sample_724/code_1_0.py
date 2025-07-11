def solve_horse_problem(n, m):
    """
    Calculates the maximum amount of water left after a horse travels m km.

    Args:
        n (int): The number of 100-liter batches of water available at the origin.
        m (float): The total distance to the destination in kilometers.
    """
    if n * 100 <= m:
        print("The horse does not have enough water to reach the destination.")
        print("<<<0>>>")
        return

    # Step 1: Find k0, the number of 100L batches at the start of the final leg.
    # This is the phase the destination 'm' falls into.
    k0 = 1  # Default k0 is 1 (final one-way trip)
    dist_boundary = 0.0
    # Iterate from phase n down to phase 2
    for k_iter in range(n, 1, -1):
        # Length of the segment to consume 100L of water
        segment_len = 100.0 / (2 * k_iter - 1)
        if m <= dist_boundary + segment_len:
            k0 = k_iter
            break
        dist_boundary += segment_len

    # Step 2: Calculate the distance covered to reach the start of the k0 phase
    # and build the string representation of the summation.
    dist_to_k0_phase = 0.0
    summation_terms_str_list = []
    # The summation runs from j = k0 + 1 up to n. If k0 = n, this loop is skipped.
    for j in range(k0 + 1, n + 1):
        term_val = 100.0 / (2 * j - 1)
        dist_to_k0_phase += term_val
        summation_terms_str_list.append(f"100/(2*{j}-1)")

    # Format the summation part of the equation string
    if not summation_terms_str_list:
        summation_str_part = "0"
    else:
        # The terms are calculated from j=n down to k0+1, but conventionally written from k0+1 to n.
        # So we reverse the list for printing to show it in ascending order of j.
        summation_str_part = f"({' + '.join(reversed(summation_terms_str_list))})"

    # Step 3: Calculate the final water amount using the derived formula.
    water_at_start_of_leg = k0 * 100.0
    dist_of_final_leg = m - dist_to_k0_phase
    consumption_rate_final_leg = 2 * k0 - 1
    water_consumed_final_leg = dist_of_final_leg * consumption_rate_final_leg
    water_left = water_at_start_of_leg - water_consumed_final_leg

    # Step 4: Construct and print the final equation.
    # This equation represents: Water_Left = Water_at_start_of_final_leg - Water_consumed_in_final_leg
    final_equation_str = (
        f"{k0}*100 - ({m} - {summation_str_part}) * (2*{k0}-1)"
    )

    print("The maximum amount of water left is given by the equation:")
    print(final_equation_str)
    print(f"= {water_left:.4f}")
    
    # Finally, print the answer in the required format.
    print(f"\n<<<{water_left:.4f}>>>")


if __name__ == '__main__':
    # --- Example Usage ---
    # n: The horse starts with n * 100 liters of water.
    # m: The distance to the destination in kilometers.
    
    # Example 1: n=3 (300L water), m=50 km
    print("--- Example 1: n=3, m=50 ---")
    solve_horse_problem(n=3, m=50.0)

    # Example 2: n=4 (400L water), m=70 km
    print("\n--- Example 2: n=4, m=70 ---")
    solve_horse_problem(n=4, m=70.0)
    
    # Example 3: n=1 (100L water), m=60km
    print("\n--- Example 3: n=1, m=60 ---")
    solve_horse_problem(n=1, m=60.0)