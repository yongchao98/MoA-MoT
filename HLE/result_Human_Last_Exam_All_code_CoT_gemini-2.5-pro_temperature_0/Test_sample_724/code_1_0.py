def solve_horse_water_problem(n, m):
    """
    Calculates the maximum amount of water left after a horse travels m km.

    Args:
        n (int): The number of 100-liter water units at the start (total n*100 liters).
        m (float): The total distance to the destination in km.
    """
    if n * 100 < m:
        print("Error: The total amount of water is less than the distance.")
        print("It's impossible to complete the journey.")
        return

    dist_traveled = 0.0

    # Iterate through the multi-trip stages, from k=n down to k=2
    for k in range(n, 1, -1):
        consumption_rate = 2 * k - 1
        # Max distance possible in this stage before water supply drops by 100L
        dist_in_stage = 100.0 / consumption_rate

        if dist_traveled + dist_in_stage >= m:
            # The destination is within this stage.
            dist_to_go = m - dist_traveled
            water_at_stage_start = k * 100.0
            water_consumed = dist_to_go * consumption_rate
            final_water = water_at_stage_start - water_consumed

            # Build the summation expression for the formula
            sum_terms = [f"1/(2*{i}-1)" for i in range(k + 1, n + 1)]
            sum_expr = " + ".join(sum_terms) if sum_terms else "0"

            print("The maximum amount of water left is given by the formula:")
            print(f"{k}*100 - ({m} - 100 * ({sum_expr})) * (2*{k}-1)")
            print(f"\nWhich evaluates to:")
            print(final_water)
            return

        # Complete the current stage and move to the next
        dist_traveled += dist_in_stage

    # If the loop completes, the destination is in the final, single-trip stage
    dist_to_go = m - dist_traveled
    water_at_stage_start = 100.0
    water_consumed = dist_to_go * 1.0  # Consumption is 1L/km
    final_water = water_at_stage_start - water_consumed

    # Build the summation expression for the formula
    sum_terms = [f"1/(2*{i}-1)" for i in range(2, n + 1)]
    sum_expr = " + ".join(sum_terms) if sum_terms else "0"

    print("The maximum amount of water left is given by the formula:")
    print(f"100 - ({m} - 100 * ({sum_expr}))")
    print(f"\nWhich evaluates to:")
    print(final_water)


if __name__ == '__main__':
    # --- User-defined inputs ---
    # n: number of 100-liter water units (e.g., n=3 means 300 liters)
    # m: distance to destination in km
    n = 3
    m = 50
    # ---------------------------

    print(f"Calculating for n={n} ({n*100} liters) and a distance m={m} km.\n")
    solve_horse_water_problem(n, m)
    # Example result for n=3, m=50:
    # Water at start of stage k=2 is 200L.
    # Distance traveled in stage k=3 is 100/(2*3-1) = 20km.
    # Remaining distance is 50-20=30km.
    # Water consumed in stage k=2 is 30km * (2*2-1)L/km = 90L.
    # Final water = 200 - 90 = 110L.
    # Formula: 2*100 - (50 - 100 * (1/(2*3-1))) * (2*2-1)
    
    # Example for n=3, m=60:
    # This will fall into the single-trip case.
    # Distance in stage k=3 is 20km. Distance in stage k=2 is 100/3 = 33.33km.
    # Total multi-trip distance = 20 + 33.33 = 53.33km.
    # Remaining distance = 60 - 53.33 = 6.67km.
    # Water at start of single-trip stage = 100L.
    # Water consumed = 6.67L.
    # Final water = 100 - 6.67 = 93.33L.
    # Formula: 100 - (60 - 100 * (1/(2*2-1) + 1/(2*3-1)))