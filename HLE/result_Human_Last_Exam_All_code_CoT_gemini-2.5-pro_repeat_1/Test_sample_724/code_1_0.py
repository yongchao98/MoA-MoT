def solve_horse_problem(n, m):
    """
    Calculates the maximum amount of water left after a horse travels m km.

    Args:
        n (int): The number of 100-liter water units available initially.
        m (float): The distance to the destination in kilometers.
    """
    if n * 100 < m:
        print("Error: The total initial water is less than the distance.")
        print("With the depot strategy, it's impossible to reach the destination.")
        # Under the simplified (but less optimal) single-trip assumption,
        # the horse can't carry all the water. If n=1, water left = 100-m.
        # But the problem implies an optimal strategy.
        return

    # Handle the simple case where only one trip is needed from the start
    if n == 1:
        if m > 100:
             print("Error: With 100L, the maximum distance is 100km.")
             return
        final_water = 100 - m
        print("With only 100L of water, the horse makes a single trip.")
        print("\nFinal Equation:")
        print(f"Remaining Water = 100 - {m}")
        print(f"Remaining Water = {final_water:.2f} L")
        return

    dist_so_far = 0.0
    summation_terms_symbolic = []
    summation_terms_numeric = []

    # Loop through phases p = 1, 2, ..., n-1
    for p in range(1, n):
        # Number of effective "loads" at the start of phase p is n-p+1
        # Consumption rate is (2 * num_loads - 1)
        rate = 2 * (n - p + 1) - 1
        
        # Max distance that can be covered in this phase before dropping another 100L
        phase_dist = 100.0 / rate

        if dist_so_far + phase_dist >= m:
            # Destination is reached in this phase (phase p)
            dist_in_this_phase = m - dist_so_far
            water_at_phase_start = (n - p + 1) * 100
            final_water = water_at_phase_start - (dist_in_this_phase * rate)
            
            print("The destination is reached before the water supply drops to 100L.")
            print("The formula for remaining water is:")
            print("(Water at phase start) - (Distance traveled in phase) * (Consumption rate)")

            sum_symbolic_str = " + ".join(summation_terms_symbolic) if summation_terms_symbolic else "0"
            sum_numeric_str = " + ".join([f"{val:.2f}" for val in summation_terms_numeric]) if summation_terms_numeric else "0"

            print("\nFinal Calculation:")
            print(f"({n}-{p}+1)*100 - ({m} - ({sum_symbolic_str})) * (2*({n}-{p}+1)-1)")
            if summation_terms_numeric:
                print(f"= {(n-p+1)*100} - ({m} - ({sum_numeric_str})) * {rate}")
                print(f"= {(n-p+1)*100} - ({m} - {dist_so_far:.2f}) * {rate}")
            else:
                print(f"= {(n-p+1)*100} - ({m} - 0) * {rate}")
            print(f"= {(n-p+1)*100} - {dist_in_this_phase:.2f} * {rate}")
            print(f"= {(n-p+1)*100} - {dist_in_this_phase * rate:.2f}")
            print(f"= {final_water:.2f} L")
            
            return

        # Phase p is completed, update distance and summation terms
        dist_so_far += phase_dist
        summation_terms_symbolic.append(f"100/(2*{n}-2*{p}+1)")
        summation_terms_numeric.append(phase_dist)

    # If the loop finishes, the destination is reached in the final leg
    dist_in_final_leg = m - dist_so_far
    final_water = 100.0 - dist_in_final_leg

    print("The horse travels until 100L of water is left, then completes the final leg.")
    print("The formula for remaining water is:")
    print("100 - (Distance traveled in final leg)")
    
    sum_symbolic_str = " + ".join(summation_terms_symbolic)
    sum_numeric_str = " + ".join([f"{val:.2f}" for val in summation_terms_numeric])

    print("\nFinal Calculation:")
    print(f"100 - ({m} - ({sum_symbolic_str}))")
    print(f"= 100 - ({m} - ({sum_numeric_str}))")
    print(f"= 100 - ({m} - {dist_so_far:.2f})")
    print(f"= 100 - {dist_in_final_leg:.2f}")
    print(f"= {final_water:.2f} L")


# --- User Input ---
# n represents the initial water supply of n * 100 liters.
# m represents the total distance in kilometers.
# Example: n=4 means 400 liters, m=100 means 100 km distance.

n = 4
m = 100

solve_horse_problem(n, m)