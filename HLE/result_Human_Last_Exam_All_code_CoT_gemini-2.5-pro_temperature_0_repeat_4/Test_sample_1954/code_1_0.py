import math

def calculate_minimax_risk(n):
    """
    Calculates the minimax risk for estimating the parameter theta of a
    Binomial(n, theta) distribution, given n i.i.d. observations.

    Args:
        n (int): The number of observations, which is also the number of trials
                 in each Binomial experiment. Must be a positive integer.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.")
        return

    # The sufficient statistic T = sum(X_i) follows a Bin(N, theta) distribution,
    # where N = n*n.
    N = n * n

    # The minimax risk for a Bin(N, theta) model under squared error loss is
    # 1 / (4 * (sqrt(N) + 1)^2).
    # Since N = n^2, sqrt(N) = n. The formula simplifies to 1 / (4 * (n + 1)^2).

    # We will show the calculation step-by-step.
    print(f"For the given n = {n}:")
    print(f"The total number of trials in the sufficient statistic is N = n * n = {N}.")
    
    # Calculate terms for the equation
    term_in_parentheses = n + 1
    squared_term = term_in_parentheses**2
    denominator = 4 * squared_term
    minimax_risk = 1 / denominator

    # Output the equation with the calculated numbers
    print("\nThe minimax risk is calculated using the formula: 1 / (4 * (n + 1)^2)")
    print("Substituting the value of n:")
    print(f"Risk = 1 / (4 * ({n} + 1)^2)")
    print(f"Risk = 1 / (4 * {term_in_parentheses}^2)")
    print(f"Risk = 1 / (4 * {squared_term})")
    print(f"Risk = 1 / {denominator}")
    print(f"\nThe final minimax risk is: {minimax_risk}")

# --- User Input ---
# You can change the value of n here.
n_value = 10
calculate_minimax_risk(n_value)