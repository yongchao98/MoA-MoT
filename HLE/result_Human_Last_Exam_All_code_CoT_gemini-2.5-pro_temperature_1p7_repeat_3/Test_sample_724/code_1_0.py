def solve_horse_water_problem():
    """
    Calculates the maximum amount of water left after a horse travels a certain distance.
    """
    # --- User-configurable values ---
    # n: The number of 100-liter containers of water at the start.
    #    Total starting water is n * 100 liters.
    n = 4
    # m: The total distance to travel in kilometers.
    m = 50

    # The problem assumes n * 100 > m, so it's possible to reach the destination.
    if n * 100 < m:
        print("Error: Not enough water to travel the distance even in a single trip.")
        return

    # This variable will track the total distance covered to reach a point
    # where the number of required trips decreases.
    dist_marker = 0.0

    # Find the segment of the journey corresponding to the distance m.
    # We loop from k=n loads down to k=2 loads.
    for k in range(n, 1, -1):
        # The length of a segment where we use up exactly 100L of water
        # while transporting 'k' loads is 100 / (2k-1) km.
        segment_len = 100.0 / (2 * k - 1)

        if m <= dist_marker + segment_len:
            # The destination 'm' is within this segment.
            # 'k' is the number of 100L loads at the start of this segment.
            
            # Build the string for the summation part of the equation.
            sum_terms_list = [f"100 / (2*{i} - 1)" for i in range(n, k, -1)]
            sum_str = " + ".join(sum_terms_list) if sum_terms_list else "0"

            # The final equation.
            water_at_start_of_segment = k * 100
            trip_factor = 2 * k - 1
            print(f"The horse has {n*100} liters of water and wants to travel {m} km.")
            print(f"The destination is reached in a segment where {k} loads are being transported.")
            print("\nThe final amount of water is given by the equation:")
            print(f"Water Left = {water_at_start_of_segment} - ({m} - ({sum_str})) * {trip_factor}")
            
            # Calculate the final result.
            water_left = k * 100.0 - (m - dist_marker) * (2 * k - 1)
            print(f"\nResult: {water_left:.4f} liters")
            return

        # Move the distance marker to the end of the completed segment.
        dist_marker += segment_len

    # If the loop finishes, the destination is in the final segment where only 1 load is left.
    # Build the string for the summation part.
    sum_terms_list = [f"100 / (2*{i} - 1)" for i in range(n, 1, -1)]
    sum_str = " + ".join(sum_terms_list)

    # The final equation for the k=1 case.
    water_at_start_of_segment = 100
    trip_factor = 1 # 2*1 - 1
    
    print(f"The horse has {n*100} liters of water and wants to travel {m} km.")
    print("The destination is reached in the final segment where only 1 load is being transported.")
    print("\nThe final amount of water is given by the equation:")
    print(f"Water Left = {water_at_start_of_segment} - ({m} - ({sum_str})) * {trip_factor}")

    # Calculate the final result.
    water_left = 100.0 - (m - dist_marker)
    print(f"\nResult: {water_left:.4f} liters")

# Execute the function to solve the problem
solve_horse_water_problem()

# Example calculation trace for n=4, m=50:
# Segment for k=4: len = 100/7 = 14.28. m > 14.28. dist_marker = 14.28
# Segment for k=3: len = 100/5 = 20.0. m (50) > 14.28+20 (34.28). dist_marker = 34.28
# Segment for k=2: len = 100/3 = 33.33. m (50) < 34.28+33.33 (67.61). We are in this segment.
# Equation will be: 2*100 - (50 - (100/(2*4-1) + 100/(2*3-1))) * (2*2-1)
# Result: 200 - (50 - (14.28 + 20)) * 3 = 200 - (15.71) * 3 = 152.857