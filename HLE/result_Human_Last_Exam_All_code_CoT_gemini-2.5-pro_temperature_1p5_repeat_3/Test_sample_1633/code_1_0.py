import math

def analyze_transformation_cost(n_values):
    """
    Analyzes the number of rewiring operations needed by calculating the
    required degree redistribution to transform a regular graph into a
    hierarchical one with the specified properties.

    The logic is as follows:
    1.  To achieve an ultra-small world with a max degree of log(n), a large number
        of "hub" nodes must be created. We estimate the number of hubs (h_n) needed
        to ensure every node is close to a hub.
    2.  We calculate the total amount of "degree" that needs to be taken from
        regular nodes and given to these hubs to raise their degree from the
        initial value (6) to the target value (log(n)).
    3.  The minimum number of rewirings, m(n), is proportional to this total
        amount of degree that must be redistributed.
    4.  We show that this value is on the order of n, meaning m(n) is in Theta(n).
    """
    print("This analysis calculates the minimum rewirings m(n) needed for the transformation.")
    print("The reasoning shows that m(n) must be on the order of n, i.e., m(n) in Theta(n).\n")

    for n in n_values:
        print(f"--- Analyzing for a network of size n = {n} ---")

        initial_degree = 6
        log_n = math.log(n)
        
        # For the model to be meaningful, n must be large enough that log(n) > 6
        if log_n <= initial_degree:
            print(f"For n={n}, log(n) is not greater than initial_degree=6. Analysis requires larger n.")
            print("="*50 + "\n")
            continue

        # Step 1: Estimate the number of hubs needed (h_n)
        # Assuming each of the (n - h_n) non-hub nodes needs to connect to one of the h_n hubs,
        # and each hub's degree is limited by log(n): h_n * log(n) >= n - h_n
        # This gives: h_n >= n / (log(n) + 1)
        num_hubs = n / (log_n + 1)

        # Step 2: Calculate total degree to be moved (Delta_degree)
        # Each hub's degree must increase from 6 to log(n)
        degree_increase_per_hub = log_n - initial_degree
        total_degree_to_move = num_hubs * degree_increase_per_hub
        
        # Step 3: m(n) is proportional to Delta_degree
        # The number of rewirings m(n) is at least the total degree to move.
        min_rewirings = total_degree_to_move
        ratio_to_n = min_rewirings / n
        
        print(f"Equation for number of hubs (h_n):")
        print(f"h_n = {n:.0f} / (log({n:.0f}) + 1) = {num_hubs:.2f}")
        print(f"Equation for total degree to move (Δ):")
        print(f"Δ = h_n * (log(n) - k₀) = {num_hubs:.2f} * ({log_n:.2f} - {initial_degree}) = {min_rewirings:.2f}")
        print("\nConclusion:")
        print(f"The minimum number of rewirings m(n) is approx. {min_rewirings:.0f}.")
        print(f"The ratio m(n) / n is approx. {ratio_to_n:.4f}.")
        print("Since this ratio approaches a constant for large n, m(n) is in Θ(n).")
        print("="*50 + "\n")


# Run the analysis for a few large values of n
analyze_transformation_cost([1000, 10000, 100000, 1000000])
<<<B>>>