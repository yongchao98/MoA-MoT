import numpy as np

def calculate_critical_correlation(Nv, Ns, mu, theta):
    """
    Calculates the critical covariance 'c_vs' required to balance potentiation and depression.

    Args:
        Nv (int): Number of neurons in input layer v.
        Ns (int): Number of neurons in input layer s.
        mu (float): Average rate of activation for v and s.
        theta (float): Heterosynaptic offset constant.

    Returns:
        float: The critical covariance c_vs.
    """
    # The formula derived from the stability analysis of the weight dynamics.
    # c_crit = mu * (mu - 1) + mu * (theta - mu) * (Nv + Ns)
    
    # Calculate the correlation coefficient rho to check if parameters are valid
    # rho must be between -1 and 1
    # rho = c_crit / (mu * (1 - mu))
    # rho = (mu * (mu - 1) + mu * (theta - mu) * (Nv + Ns)) / (mu * (1 - mu))
    # rho = -1 + (theta - mu) * (Nv + Ns) / (1 - mu)
    
    rho = -1 + (theta - mu) * (Nv + Ns) / (1 - mu)
    
    if not (-1 <= rho <= 1):
        print("Warning: The given parameters result in a correlation coefficient outside the valid range [-1, 1].")
        print(f"Calculated correlation coefficient (rho): {rho:.4f}")
        print("This means a balance point with this level of correlation is not physically possible.")
        print("Adjust parameters, e.g., ensure (theta - mu) * (Nv + Ns) is between 0 and 2 * (1 - mu).\n")


    c_crit = mu * (mu - 1) + mu * (theta - mu) * (Nv + Ns)
    return c_crit

if __name__ == '__main__':
    # Example parameters
    # Nv and Ns are the number of neurons in the respective input layers.
    Nv = 100
    # Ns is the number of neurons in the s layer
    Ns = 100
    # mu is the average firing rate, which should be a small value.
    mu = 0.01
    # theta is the depression threshold. For balance, theta should be slightly larger than mu.
    theta = 0.015

    # Calculate the critical correlation
    critical_covariance = calculate_critical_correlation(Nv, Ns, mu, theta)

    # Output the result, showing the numbers used in the final equation as requested.
    print("The formula for the critical covariance (c_crit) is:")
    print("c_crit = mu*(mu - 1) + mu*(theta - mu)*(Nv + Ns)\n")
    
    print("Plugging in the example values:")
    print(f"c_crit = {mu}*({mu} - 1) + {mu}*({theta} - {mu})*({Nv} + {Ns})")
    
    # Calculate intermediate terms for clarity
    term1 = mu * (mu - 1)
    term2 = mu * (theta - mu) * (Nv + Ns)
    print(f"c_crit = {term1:.4f} + {term2:.4f}")

    print("\nFinal Result:")
    print(f"The critical amount of correlation (covariance) is: {critical_covariance:.4f}")
