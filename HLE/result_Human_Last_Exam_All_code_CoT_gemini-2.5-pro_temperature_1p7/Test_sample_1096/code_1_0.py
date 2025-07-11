import math

def calculate_optimal_steps(M):
    """
    Calculates the optimal two-step gradient descent learning rates.

    Args:
        M (float): The conditioning number kappa of the function.
                   We assume m=1, so M=kappa. Must be >= 1.

    Returns:
        tuple: A tuple containing the two optimal learning rates (gamma1, gamma2).
    """
    if M < 1:
        print("Error: M (kappa) must be greater than or equal to 1.")
        return None, None
    
    print(f"Given M = kappa = {M}")
    
    # Using the formulas derived from the Chebyshev polynomial optimization:
    # gamma_1,2 = (4*(M+1) -+ 2*sqrt(2)*(M-1)) / (M^2 + 6*M + 1)
    
    # Calculate the denominator
    denominator = M**2 + 6*M + 1
    
    # Calculate the two terms in the numerator
    term1_num = 4 * (M + 1)
    term2_num = 2 * math.sqrt(2) * (M - 1)
    
    # Calculate the full numerators for gamma1 and gamma2
    numerator1 = term1_num - term2_num
    numerator2 = term1_num + term2_num
    
    # Calculate gamma1 and gamma2
    gamma1 = numerator1 / denominator
    gamma2 = numerator2 / denominator
    
    print("\nThe optimal step sizes are calculated as follows:")
    print("Formula: (4*(M+1) \u2213 2*sqrt(2)*(M-1)) / (M\u00b2 + 6*M + 1)\n")

    # Print the calculation for gamma1
    print("For gamma_1 (with the minus sign):")
    print(f"  Numerator   = 4*({M}+1) - 2*sqrt(2)*({M}-1) = {term1_num:.4f} - {term2_num:.4f} = {numerator1:.4f}")
    print(f"  Denominator = {M}\u00b2 + 6*{M} + 1 = {denominator:.4f}")
    print(f"  gamma_1     = {numerator1:.4f} / {denominator:.4f} = {gamma1:.4f}\n")

    # Print the calculation for gamma2
    print("For gamma_2 (with the plus sign):")
    print(f"  Numerator   = 4*({M}+1) + 2*sqrt(2)*({M}-1) = {term1_num:.4f} + {term2_num:.4f} = {numerator2:.4f}")
    print(f"  Denominator = {M}\u00b2 + 6*{M} + 1 = {denominator:.4f}")
    print(f"  gamma_2     = {numerator2:.4f} / {denominator:.4f} = {gamma2:.4f}\n")

    print(f"The best choice for the pair (gamma_1, gamma_2) is ({gamma1:.4f}, {gamma2:.4f}).")
    
    return gamma1, gamma2

if __name__ == '__main__':
    # You can change this value to any conditioning number M >= 1
    M_kappa = 10.0
    calculate_optimal_steps(M_kappa)