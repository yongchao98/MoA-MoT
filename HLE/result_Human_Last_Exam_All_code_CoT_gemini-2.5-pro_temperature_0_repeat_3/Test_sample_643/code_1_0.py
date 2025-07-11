def calculate_wet_trip_fraction(n, p):
    """
    Calculates the asymptotic fraction of trips where the professor gets wet.

    Args:
        n (int): The total number of umbrellas.
        p (float): The probability of rain for any given trip, between 0 and 1.

    Returns:
        float: The fraction of trips the professor gets wet.
    """
    if not (isinstance(n, int) and n > 0):
        print("Error: n must be a positive integer.")
        return
    if not (isinstance(p, float) and 0 < p < 1):
        print("Error: p must be a float between 0 and 1.")
        return

    numerator = p * (1 - p)
    denominator = n + 1 - p
    fraction = numerator / denominator

    print(f"Given n = {n} umbrellas and a probability of rain p = {p}:")
    print("The formula for the fraction of wet trips is: p * (1-p) / (n + 1 - p)")
    print(f"Plugging in the values:")
    print(f"= {p} * (1 - {p}) / ({n} + 1 - {p})")
    print(f"= {p} * {1-p} / {n + 1 - p}")
    print(f"= {numerator} / {denominator}")
    print(f"= {fraction}")

# Example usage with n=5 umbrellas and p=0.2 probability of rain.
calculate_wet_trip_fraction(n=5, p=0.2)
