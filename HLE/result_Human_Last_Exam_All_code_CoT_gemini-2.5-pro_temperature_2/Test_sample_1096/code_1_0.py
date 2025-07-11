import math

def calculate_optimal_steps(M):
    """
    Calculates the optimal step sizes (gamma_1, gamma_2) for two-step gradient descent.

    Args:
        M (float): The condition number (kappa) of the function. Must be > 1.

    Returns:
        tuple: A tuple containing the two optimal step sizes (gamma_1, gamma_2).
    """
    if M <= 1:
        raise ValueError("M (condition number) must be greater than 1.")

    # The formulas are derived from Chebyshev polynomial theory.
    # The optimal step sizes are the roots of the quadratic equation z^2 - s*z + p = 0,
    # where s and p are the sum and product of the steps.
    # s = (gamma_1 + gamma_2) = 8*(M+1) / (M^2 + 6*M + 1)
    # p = (gamma_1 * gamma_2) = 8 / (M^2 + 6*M + 1)
    
    # Direct formulas for the roots:
    # gamma_{1,2} = (4*(M+1) +/- 2*sqrt(2)*(M-1)) / (M**2 + 6*M + 1)
    
    denominator = M**2 + 6*M + 1
    term1 = 4 * (M + 1)
    term2 = 2 * math.sqrt(2) * (M - 1)
    
    gamma1 = (term1 - term2) / denominator
    gamma2 = (term1 + term2) / denominator
    
    return gamma1, gamma2

if __name__ == '__main__':
    # We are given M = kappa, m=1
    # Example usage:
    try:
        M_input_str = input("Enter the condition number M (kappa): ")
        M = float(M_input_str)
        
        gamma1, gamma2 = calculate_optimal_steps(M)

        print(f"For a condition number M = {M}:")
        print(f"The best choice for the pair (gamma_1, gamma_2) is:")
        print(f"gamma_1 = {gamma1}")
        print(f"gamma_2 = {gamma2}")

    except ValueError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
