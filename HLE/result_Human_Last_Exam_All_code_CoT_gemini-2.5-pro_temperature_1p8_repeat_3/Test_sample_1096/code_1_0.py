import math

def compute_optimal_gammas(M):
    """
    Computes the optimal 2-step gradient descent learning rates (gammas).

    Args:
        M (float): The condition number kappa, which is the ratio of the smoothness
                   constant M to the strong convexity constant mu (assumed to be 1).

    Returns:
        tuple: A tuple containing the optimal (gamma_1, gamma_2).
    """
    if M < 1:
        raise ValueError("Condition number M must be >= 1.")
    if M == 1:
        # For M=1, the function is perfectly conditioned, problem is trivial.
        # The formula would lead to division by zero, but any reasonable step works.
        # The optimal step is 1/M = 1.
        return (1.0, 1.0)

    # Based on the Chebyshev polynomial method, the inverses of the optimal
    # gammas are the roots of the optimal polynomial.
    # The roots are given by the formula: (M+1)/2 +/- (M-1)/(2*sqrt(2))
    val_sqrt2 = math.sqrt(2.0)
    
    # Let's calculate the two roots which are 1/gamma_1 and 1/gamma_2
    root1 = (M + 1.0) / 2.0 - (M - 1.0) / (2.0 * val_sqrt2)
    root2 = (M + 1.0) / 2.0 + (M - 1.0) / (2.0 * val_sqrt2)
    
    # The optimal gammas are the inverses of these roots.
    # The order is interchangeable.
    gamma1 = 1.0 / root1
    gamma2 = 1.0 / root2
    
    return gamma1, gamma2

# Set the condition number M (kappa). Let's use M=10 as an example.
M = 10.0

# Calculate the optimal gamma pair
try:
    gamma1, gamma2 = compute_optimal_gammas(M)

    # --- Outputting the results as requested ---
    print(f"For a condition number M = {M}:")
    
    # Explaining the formula for gamma_1
    print("\nThe best choice for gamma_1 is derived from one of the roots of the optimal polynomial:")
    print(f"1/gamma_1 = ({M} + 1)/2 - ({M} - 1)/(2 * sqrt(2))")
    print(f"gamma_1 = {gamma1}")

    # Explaining the formula for gamma_2
    print("\nThe best choice for gamma_2 is derived from the other root:")
    print(f"1/gamma_2 = ({M} + 1)/2 + ({M} - 1)/(2 * sqrt(2))")
    print(f"gamma_2 = {gamma2}")

except ValueError as e:
    print(e)
