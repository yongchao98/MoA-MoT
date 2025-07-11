import math

def calculate_optimal_stepsizes(M):
    """
    Calculates the optimal two-step gradient descent learning rates (gamma_1, gamma_2)
    for M-smooth and 1-strongly convex functions.

    Args:
        M (float): The condition number (M/mu), where mu is assumed to be 1.
    """
    if M <= 1:
        print("M (the condition number) must be greater than 1.")
        return

    # Denominator D = M^2 + 6M + 1
    D_num = M**2 + 6*M + 1
    
    # Numerators for the two gammas
    # N = 4(M+1) +- 2*sqrt(2)*(M-1)
    term1_num = 4 * (M + 1)
    term2_num = 2 * math.sqrt(2) * (M - 1)
    
    # Calculate gamma_1 and gamma_2
    gamma1_num = term1_num - term2_num
    gamma2_num = term1_num + term2_num
    
    gamma1 = gamma1_num / D_num
    gamma2 = gamma2_num / D_num
    
    print(f"For M = {M}:")
    print("-" * 20)
    print("The optimal step sizes (gamma_1, gamma_2) are given by the formulas:")
    print("gamma_1 = (4 * (M + 1) - 2 * sqrt(2) * (M - 1)) / (M^2 + 6*M + 1)")
    print("gamma_2 = (4 * (M + 1) + 2 * sqrt(2) * (M - 1)) / (M^2 + 6*M + 1)")
    print("\n" + "-"*20)

    # Printing the breakdown of the calculation for gamma_1
    print(f"For gamma_1:")
    print(f"  Numerator   = 4 * ({M} + 1) - 2 * sqrt(2) * ({M} - 1) = {term1_num:.4f} - {term2_num:.4f} = {gamma1_num:.4f}")
    print(f"  Denominator = {M}^2 + 6*{M} + 1 = {M**2} + {6*M} + 1 = {D_num:.4f}")
    print(f"  gamma_1 = {gamma1_num:.4f} / {D_num:.4f} = {gamma1:.4f}\n")
    
    # Printing the breakdown of the calculation for gamma_2
    print(f"For gamma_2:")
    print(f"  Numerator   = 4 * ({M} + 1) + 2 * sqrt(2) * ({M} - 1) = {term1_num:.4f} + {term2_num:.4f} = {gamma2_num:.4f}")
    print(f"  Denominator = {M}^2 + 6*{M} + 1 = {M**2} + {6*M} + 1 = {D_num:.4f}")
    print(f"  gamma_2 = {gamma2_num:.4f} / {D_num:.4f} = {gamma2:.4f}")


# Example usage with M = 10
M_example = 10
calculate_optimal_stepsizes(M_example)