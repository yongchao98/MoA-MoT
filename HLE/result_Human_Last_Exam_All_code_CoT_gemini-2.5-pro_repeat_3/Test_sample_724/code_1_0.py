def solve_water_problem(n, m):
    """
    Calculates the maximum amount of water left after a horse travels m km.
    This function implements the optimal caching strategy for the "Jeep problem".

    Args:
        n (int): The number of 100L water loads available at the origin.
        m (float): The total distance to the destination in km.
    """
    # --- Input Validation ---
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.")
        return
    if not isinstance(m, (int, float)) or m < 0:
        print("Error: m must be a non-negative number.")
        return
    if n * 100 <= m:
        # This is a simplification. The actual max distance is a harmonic-like series.
        # But if total water is less than the distance, it's clearly impossible.
        print(f"Error: Not enough water to reach the destination.")
        print(f"Total water {n*100}L is not sufficient to travel {m}km.")
        return

    # --- Calculation ---
    d_cumulative = 0.0
    sum_terms = []

    # Loop through the n-1 caching legs
    for k in range(1, n):
        num_loads = n - k + 1
        rate = float(2 * num_loads - 1)
        leg_distance = 100.0 / rate

        if m <= d_cumulative + leg_distance:
            # Destination is in this intermediate leg
            distance_in_leg = m - d_cumulative
            water_left = (num_loads * 100.0) - (distance_in_leg * rate)

            # Build the equation string for the amount of water left
            if not sum_terms:
                d_prev_str = "0"
            else:
                d_prev_str = f"100 * ({' + '.join(sum_terms)})"

            equation = f"{num_loads} * 100 - ({m} - {d_prev_str}) * {int(rate)}"

            print("Final Equation (in terms of summation):")
            print(f"Water Left = {equation}")

            print("\nCalculation:")
            print(f"Water Left = {num_loads * 100} - ({m} - {d_cumulative:.2f}) * {int(rate)}")
            print(f"             = {num_loads * 100} - {distance_in_leg:.2f} * {int(rate)}")
            print(f"             = {num_loads * 100} - {distance_in_leg * rate:.2f}")
            print(f"             = {water_left:.2f}")

            print(f"\nThe maximum amount of water left is {water_left:.2f} liters.")
            print(f"<<<{water_left:.2f}>>>")
            return

        d_cumulative += leg_distance
        sum_terms.append(f"1 / {int(rate)}")

    # If the loop completes, the destination is in the final leg (with 1 load)
    distance_in_final_leg = m - d_cumulative
    water_left = 100.0 - distance_in_final_leg

    d_prev_str = f"100 * ({' + '.join(sum_terms)})"
    equation = f"100 - ({m} - {d_prev_str})"

    print("Final Equation (in terms of summation):")
    print(f"Water Left = {equation}")

    print("\nCalculation:")
    print(f"Water Left = 100 - ({m} - {d_cumulative:.2f})")
    print(f"             = 100 - {distance_in_final_leg:.2f}")
    print(f"             = {water_left:.2f}")

    print(f"\nThe maximum amount of water left is {water_left:.2f} liters.")
    print(f"<<<{water_left:.2f}>>>")


if __name__ == '__main__':
    # --- Set your inputs here ---
    # n: number of 100-liter water loads
    # m: distance to travel in kilometers
    n = 3
    m = 60

    solve_water_problem(n, m)