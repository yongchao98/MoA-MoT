def calculate_expected_ratio_limit(n_max):
    """
    Calculates the expected number of remaining numbers, E_n,
    and the ratio E_n/n up to n_max.
    """
    if n_max < 2:
        print("Please provide n_max >= 2.")
        return

    # E[n] will store the expected number of remaining items for size n.
    E = [0.0] * (n_max + 1)
    
    # Base cases
    # E[0] = 0 (already initialized)
    # For a single number, it remains.
    if n_max >= 1:
        E[1] = 1.0
    # For two numbers (1, 2), they are crossed out.
    # E[2] = 0 (already initialized)

    # Use the recurrence relation to compute E[n] for n >= 3
    # E[n] = ((n-2)/(n-1)) * E[n-1] + (2/(n-1)) * E[n-2]
    for n in range(3, n_max + 1):
        term1_numerator = n - 2
        term1_denominator = n - 1
        term2_numerator = 2
        term2_denominator = n - 1
        
        E[n] = (term1_numerator / term1_denominator) * E[n-1] + (term2_numerator / term2_denominator) * E[n-2]

    # The final equation is Ratio = E[n] / n
    final_n = n_max
    final_En = E[n_max]
    final_ratio = final_En / final_n
    
    print(f"For n = {final_n}:")
    print(f"The expected number of remaining numbers, E_n, is: {final_En}")
    print("The final equation for the ratio is:")
    print(f"Ratio = {final_En} / {final_n}")
    print(f"The value of the ratio is approximately: {final_ratio}")
    print("\nAs n approaches infinity, this ratio converges to 1/3.")

# Set a large value for n to approximate the limit
n_large = 10000
calculate_expected_ratio_limit(n_large)