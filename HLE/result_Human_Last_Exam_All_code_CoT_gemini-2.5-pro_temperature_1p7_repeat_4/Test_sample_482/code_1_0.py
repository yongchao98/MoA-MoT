import math

def calculate_critical_correlation(mu, theta, N):
    """
    Calculates the critical correlation C = E[v s] required to balance
    potentiation and depression in the network, based on the derived formula.
    
    The function also prints the formula and a step-by-step calculation
    with the provided numerical values, as requested.

    Args:
        mu (float): The average rate of activation for input neurons (μ).
        theta (float): The heterosynaptic offset constant (θ).
        N (int): The number of neurons in each input population (assuming N_v = N_s = N).

    Returns:
        float: The value of the critical correlation C.
    """
    if not (isinstance(mu, (int, float)) and isinstance(theta, (int, float)) and isinstance(N, int)):
        raise TypeError("Inputs mu and theta must be numeric, and N must be an integer.")
    if not (0 <= mu <= 1):
        print("Warning: mu (as a rate from a Binomial process) should typically be between 0 and 1.")
    if N <= 0:
        raise ValueError("N must be a positive integer.")

    # Derived formula for the critical correlation C = E[v*s]
    critical_C = 2 * mu * theta - mu**2 - (mu - mu**2) / N

    print("The derived formula for the critical correlation, C = E[v s], is:")
    # Using Greek letters for a clear representation of the formula
    print(f"C = 2 * \u03BC * \u03B8 - \u03BC\u00B2 - (\u03BC - \u03BC\u00B2) / N\n")
    
    print("--- Calculation with Provided Values ---")
    print(f"Given values are:")
    print(f"Average rate \u03BC = {mu}")
    print(f"Offset constant \u03B8 = {theta}")
    print(f"Population size N = {N}\n")

    # Displaying the calculation step-by-step with the numbers
    print("Substituting the values into the equation:")
    print(f"C = 2 * {mu} * {theta} - {mu}**2 - ({mu} - {mu}**2) / {N}")
    
    # Calculate each part of the equation
    term1 = 2 * mu * theta
    term2 = mu**2
    numerator = mu - mu**2
    term3 = numerator / N

    print(f"C = {term1} - {term2} - ({numerator}) / {N}")
    print(f"C = {term1 - term2} - {term3}")
    
    print("\nFinal Result:")
    print(f"The critical amount of correlation C = {critical_C}")
    
    return critical_C

# Since no specific numerical values were given in the problem description,
# we use placeholder values to demonstrate the function's usage.
# You can change these values to explore different scenarios.
example_mu = 0.1
example_theta = 0.2
example_N = 100

# Execute the calculation
calculate_critical_correlation(mu=example_mu, theta=example_theta, N=example_N)
