import math

def calculate_minimax_risk(n):
    """
    Calculates the minimax risk for estimating the parameter theta of a Binomial distribution.

    Args:
        n (int): The number of i.i.d. observations from a Bin(n, theta) distribution.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.")
        return

    # Based on the problem statement, we have n i.i.d. observations X_i ~ Bin(n, theta).
    # The sufficient statistic is T = sum(X_i), which follows Bin(n*n, theta).
    # Let N be the total number of trials.
    N = n * n

    # The general formula for the minimax risk for Bin(N, theta) with squared error loss
    # is 1 / (4 * (1 + sqrt(N))^2).
    # Substituting N = n^2, the formula simplifies to 1 / (4 * (1 + n)^2).

    denominator_val = 4 * (1 + n)**2
    risk_value = 1 / denominator_val

    print(f"For n = {n}:")
    print(f"The sufficient statistic is the sum of observations, which follows a Bin(N, theta) distribution with N = n*n = {N}.")
    print("The minimax risk is calculated using the formula: 1 / (4 * (1 + n)^2)")
    print(f"The equation with the given n is: 1 / (4 * (1 + {n})^2) = 1 / {denominator_val}")
    print(f"The numerical value of the minimax risk is: {risk_value}")


# Example usage with n = 10
calculate_minimax_risk(10)