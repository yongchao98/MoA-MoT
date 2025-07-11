import numpy as np

def calculate_critical_correlation(mu, theta):
    """
    Calculates the critical amount of correlation (covariance) required
    to balance potentiation and depression in the described neural network.

    Args:
        mu (float): The average rate of activation for input populations.
        theta (float): The heterosynaptic offset constant.

    Returns:
        float: The critical covariance c_crit.
    """
    if theta <= mu:
        print("Warning: For selectivity to occur, theta should be greater than mu.")
        print("In this case, the critical covariance is zero or negative, implying any positive correlation leads to non-selective growth.")
        return 0

    # The formula for the critical covariance is c_crit = mu * (theta - mu)
    c_crit = mu * (theta - mu)
    
    # Output the explanation and the final equation with values
    print("The critical amount of correlation (covariance) is determined by the balance between Hebbian potentiation and heterosynaptic depression.")
    print("This balance is reached when the covariance 'c' satisfies the equation: c + mu^2 = theta * mu")
    print("Solving for the critical covariance, c_crit, we get the formula: c_crit = mu * (theta - mu)")
    print("\nUsing the provided example values:")
    print(f"Average activation rate (mu) = {mu}")
    print(f"Heterosynaptic offset (theta) = {theta}")
    print("\nThe calculation is:")
    print(f"c_crit = {mu} * ({theta} - {mu})")
    print(f"c_crit = {mu} * ({theta - mu})")
    print(f"c_crit = {c_crit}")
    
    return c_crit

# Example values for the parameters
# The average rate of activation (e.g., in spikes/sec)
mu = 5.0
# The heterosynaptic offset constant (must be > mu for this scenario)
theta = 8.0

# Calculate and print the result
critical_covariance = calculate_critical_correlation(mu, theta)