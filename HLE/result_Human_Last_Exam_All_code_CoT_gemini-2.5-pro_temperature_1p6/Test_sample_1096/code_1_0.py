import math

def find_optimal_steps(M):
    """
    Calculates the optimal two-step gradient descent learning rates (gammas)
    for an M-smooth and 1-strongly convex function.

    Args:
        M (float): The condition number (smoothness parameter / strong convexity parameter).

    Returns:
        tuple: A pair of optimal learning rates (gamma_1, gamma_2).
    """
    if M < 1:
        raise ValueError("Condition number M must be >= 1.")
    if M == 1:
        # For M=1, the function is quadratic with a single eigenvalue.
        # One step is enough to reach the minimum with gamma = 1/M = 1.
        # For two steps, the problem is ill-defined, but we can return the limit.
        return (1.0, 1.0)
        
    kappa = M

    # The optimal step sizes are derived from the roots of the shifted Chebyshev polynomial of degree 2.
    # The formulas for the step sizes gamma_1 and gamma_2 are:
    # gamma_{1,2} = (4*(kappa+1) +- 2*sqrt(2)*(kappa-1)) / (kappa^2 + 6*kappa + 1)
    
    # Numerator parts
    num_term1 = 4 * (kappa + 1)
    num_term2 = 2 * math.sqrt(2) * (kappa - 1)
    
    # Denominator
    den = kappa**2 + 6*kappa + 1
    
    # Calculate gamma_1 and gamma_2
    gamma1_num = num_term1 + num_term2
    gamma1 = gamma1_num / den
    
    gamma2_num = num_term1 - num_term2
    gamma2 = gamma2_num / den

    # Output the explanation and the final formulas with numerical values.
    print(f"For a condition number M = {M}:")
    print("-" * 30)
    
    print("The best choice for the pair (gamma_1, gamma_2) is given by the formulas:")
    print("gamma_{1,2} = (4*(M+1) \u00B1 2*sqrt(2)*(M-1)) / (M^2 + 6*M + 1)\n")

    print("gamma_1:")
    print(f"  = (4*({kappa}+1) + 2*sqrt(2)*({kappa}-1)) / ({kappa}^2 + 6*{kappa} + 1)")
    print(f"  = ({num_term1:.4f} + {num_term2:.4f}) / ({den:.4f})")
    print(f"  = {gamma1_num:.4f} / {den:.4f}")
    print(f"  = {gamma1:.6f}\n")

    print("gamma_2:")
    print(f"  = (4*({kappa}+1) - 2*sqrt(2)*({kappa}-1)) / ({kappa}^2 + 6*{kappa} + 1)")
    print(f"  = ({num_term1:.4f} - {num_term2:.4f}) / ({den:.4f})")
    print(f"  = {gamma2_num:.4f} / {den:.4f}")
    print(f"  = {gamma2:.6f}\n")
    
    print("Note: The provided definition S=sqrt(M^2+(M-1)^2) is not used in this standard derivation and might be a distractor or belong to a different context.")

    return gamma1, gamma2

# Example usage with M = 10 (a common choice for demonstration)
# The condition number kappa = M
M = 10
find_optimal_steps(M)
