def solve_horse_problem(n, m):
    """
    Calculates the maximum amount of water left after a horse trip.

    Args:
        n (int): The number of 100-liter units of water at the start.
        m (float): The total distance to travel in kilometers.
    """

    # --- Step 1: Find the final phase number K ---
    # A journey ends in phase K. We find K by determining where the total
    # distance 'm' falls within the cumulative phase distances.
    dist_covered_before_phase = 0.0
    final_phase_K = 1  # Default to 1 for the last leg
    for k in range(n, 1, -1):
        # The max distance that can be covered in phase k is 100L / (2k-1)L/km
        phase_dist = 100.0 / (2 * k - 1)
        if m <= dist_covered_before_phase + phase_dist:
            final_phase_K = k
            break
        dist_covered_before_phase += phase_dist
    
    K = final_phase_K

    # --- Step 2: Build the summation string and calculate its value ---
    # The summation represents the total distance covered in phases > K.
    summation_terms = []
    summation_value = 0.0
    # The loop goes from n down to K+1
    for i in range(n, K, -1):
        term_str = f"100/(2*{i}-1)"
        summation_terms.append(term_str)
        summation_value += 100.0 / (2 * i - 1)
    
    # Create a string representation for the summation part of the equation
    summation_str = " + ".join(summation_terms)
    if not summation_str:
        summation_str = "0"

    # --- Step 3: Construct the final equation and calculate the result ---
    # The formula for water left is: K*100 - (2K-1) * (m - sum_of_distances)
    final_equation = f"{K}*100 - (2*{K}-1) * ({m} - ({summation_str}))"
    water_left = K * 100.0 - (2 * K - 1) * (m - summation_value)

    # --- Step 4: Print the results ---
    print(f"For n={n} (initial water = {n*100}L) and destination at m={m}km:")
    print(f"The journey will end during phase K={K}.")
    print("\nThe general formula for the maximum water remaining is:")
    print("Water Left = K*100 - (2*K - 1) * (m - (Sum_{i=K+1 to n} [100/(2i-1)]))")

    print("\nPlugging the numbers into the formula, we get the specific equation:")
    print(final_equation)
    
    print(f"\nCalculated result:")
    print(f"Maximum water left at destination = {water_left:.4f} liters")


# --- User-configurable inputs ---
# n: Initial water is n * 100 liters
# m: Distance to the destination in km
n = 4
m = 100

# Solve the problem with the given inputs
solve_horse_problem(n, m)
