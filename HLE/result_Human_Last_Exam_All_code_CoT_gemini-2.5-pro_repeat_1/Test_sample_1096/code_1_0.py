import math

def calculate_optimal_step_sizes(M, mu=1):
    """
    Calculates the optimal two-step gradient descent learning rates.

    Args:
        M (float): The smoothness constant (often denoted as kappa, the condition number, when mu=1).
        mu (float): The strong convexity constant.

    Returns:
        tuple: A pair of optimal learning rates (gamma_1, gamma_2).
    """
    if M < mu:
        raise ValueError("M must be greater than or equal to mu.")
    if M == mu:
        # If M=mu, the problem is perfectly conditioned.
        # Converges in one step with gamma = 1/M.
        return (1/M, 1/M)

    # The optimal step sizes are the roots of a specific quadratic equation
    # derived from a shifted Chebyshev polynomial. The formulas are:
    # gamma_1 = 2 / ((M + mu) + (M - mu) / sqrt(2))
    # gamma_2 = 2 / ((M + mu) - (M - mu) / sqrt(2))
    
    sqrt2 = math.sqrt(2)
    
    denominator1 = (M + mu) + (M - mu) / sqrt2
    gamma1 = 2 / denominator1
    
    denominator2 = (M + mu) - (M - mu) / sqrt2
    gamma2 = 2 / denominator2
    
    return (gamma1, gamma2)

if __name__ == '__main__':
    # Per the problem statement, we assume mu=1 and M=kappa.
    # Let's use a sample value for the condition number, e.g., M = 10.
    M_kappa = 10
    mu = 1

    print(f"Calculating optimal step sizes for M = {M_kappa} and mu = {mu}:")
    
    try:
        gamma1, gamma2 = calculate_optimal_step_sizes(M_kappa, mu)
        
        # The final output should show the equation with the numbers.
        print("\nThe optimal pair of step sizes (gamma_1, gamma_2) is:")
        print(f"gamma_1 = {gamma1}")
        print(f"gamma_2 = {gamma2}")
        print("\nThe order of the step sizes does not matter.")

    except ValueError as e:
        print(f"Error: {e}")
