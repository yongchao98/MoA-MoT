import math

def solve_optimal_steps(M):
    """
    Calculates the optimal two-step gradient descent step sizes (gamma_1, gamma_2)
    for an M-smooth and 1-strongly convex function.
    
    Args:
        M (float): The condition number (kappa).
    """
    if M < 1:
        print("M must be >= 1.")
        return
    if M == 1:
        print("The function is quadratic. Optimal step is 1. One step is enough.")
        gamma1 = 1.0
        gamma2 = 0.0 # Not needed
        print(f"For M = {M}:")
        print(f"gamma_1 = {gamma1:.6f}")
        print(f"gamma_2 = {gamma2:.6f}")
        return

    # Denominator for the coefficients
    denominator = M**2 + 6*M + 1

    # Sum of the gammas
    sum_gammas = (8 * (M + 1)) / denominator

    # Product of the gammas
    prod_gammas = 8 / denominator
    
    print(f"For M = {M}:")
    print("The optimal step sizes gamma_1 and gamma_2 are the roots of the quadratic equation:")
    print(f"x^2 - ({sum_gammas:.6f}) * x + ({prod_gammas:.6f}) = 0")
    print("\nThe final equation is:")
    print(f"x^2 - (8 * (M+1) / (M^2+6*M+1)) * x + (8 / (M^2+6*M+1)) = 0")
    print("\nWhere each number in the equation is:")
    print(f"Coefficient of x^2: 1")
    print(f"Coefficient of x: -{sum_gammas:.6f}")
    print(f"Constant term: {prod_gammas:.6f}\n")


    # Solving the quadratic equation x^2 - sum_gammas*x + prod_gammas = 0
    # Discriminant
    discriminant = sum_gammas**2 - 4 * prod_gammas

    # Roots
    sqrt_discriminant = math.sqrt(discriminant)
    gamma1 = (sum_gammas - sqrt_discriminant) / 2
    gamma2 = (sum_gammas + sqrt_discriminant) / 2
    
    print("The optimal step sizes are:")
    print(f"gamma_1 = {gamma1:.6f}")
    print(f"gamma_2 = {gamma2:.6f}")


# Example usage with M = 10 (as kappa)
M = 10
solve_optimal_steps(M)
