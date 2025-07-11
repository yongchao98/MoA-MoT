import math

def find_optimal_gammas(M):
    """
    Calculates the optimal two-step gradient descent learning rates (gammas)
    for an M-smooth and 1-strongly convex function.

    Args:
        M (float): The smoothness constant (conditioning number, kappa), M > 1.

    Returns:
        tuple: A pair of floats (gamma1, gamma2) representing the optimal learning rates.
    """
    if M <= 1:
        raise ValueError("M must be greater than 1.")

    # The optimal gamma values are the roots of the quadratic equation:
    # (M^2 + 6M + 1) * z^2 - 8(M+1) * z + 8 = 0
    # We solve this using the quadratic formula.
    
    # Denominator of the gamma expressions
    denominator = M**2 + 6*M + 1
    
    # Common term in the numerator
    term1_num = 4 * (M + 1)
    
    # Second term in the numerator, involving the square root
    term2_num = 2 * math.sqrt(2) * (M - 1)
    
    # Calculate the two roots for gamma
    gamma1 = (term1_num + term2_num) / denominator
    gamma2 = (term1_num - term2_num) / denominator
    
    return gamma1, gamma2

def main():
    """
    Main function to demonstrate the calculation for a sample value of M.
    """
    # We use M=kappa as per the problem description.
    # Let's use an example value for M, e.g., M=10.
    M = 10.0
    
    print(f"For a function with smoothness M = {M} and strong convexity mu = 1:")
    
    gamma1, gamma2 = find_optimal_gammas(M)
    
    # The final equations for gamma1 and gamma2 are:
    # gamma = (4*(M+1) +/- 2*sqrt(2)*(M-1)) / (M^2 + 6*M + 1)
    
    denominator = M**2 + 6*M + 1
    term1_num = 4 * (M + 1)
    term2_num = 2 * math.sqrt(2) * (M - 1)

    print("\nThe best choice for the pair (gamma_1, gamma_2) is given by the formulas:")
    print("gamma_1 = (4*(M+1) + 2*sqrt(2)*(M-1)) / (M^2 + 6*M + 1)")
    print("gamma_2 = (4*(M+1) - 2*sqrt(2)*(M-1)) / (M^2 + 6*M + 1)")
    print("\nFor M = 10:")
    print(f"Denominator = {M**2} + 6*{M} + 1 = {denominator}")
    print(f"Numerator term 1 = 4*({M}+1) = {term1_num}")
    print(f"Numerator term 2 = 2*sqrt(2)*({M}-1) = {term2_num:.4f}")
    
    print("\nThe numerical values are:")
    print(f"gamma_1 = ({term1_num:.4f} + {term2_num:.4f}) / {denominator:.4f} = {gamma1:.4f}")
    print(f"gamma_2 = ({term1_num:.4f} - {term2_num:.4f}) / {denominator:.4f} = {gamma2:.4f}")


if __name__ == "__main__":
    main()
