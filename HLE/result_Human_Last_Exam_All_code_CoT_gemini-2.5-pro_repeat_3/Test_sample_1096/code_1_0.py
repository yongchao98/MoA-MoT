import math

def find_optimal_steps(M):
    """
    Calculates the optimal two-step gradient descent learning rates.

    Args:
        M (float): The condition number (kappa), M >= 1.

    Returns:
        tuple: A pair of floats (gamma_1, gamma_2) representing the optimal learning rates.
    """
    if M < 1:
        raise ValueError("Condition number M must be >= 1")
    if M == 1:
        # For M=1, the problem is trivial and a single step with gamma=1 suffices.
        # The formula gives gamma1=gamma2=1.0
        return (1.0, 1.0)

    # Denominator of the expression for gamma_1 and gamma_2
    denominator = M**2 + 6*M + 1

    # The term under the square root in the quadratic formula
    sqrt_term = 2 * math.sqrt(2) * (M - 1)

    # Numerator parts
    common_numerator_part = 4 * (M + 1)

    # Calculate gamma_1 and gamma_2
    gamma_1 = (common_numerator_part - sqrt_term) / denominator
    gamma_2 = (common_numerator_part + sqrt_term) / denominator
    
    return gamma_1, gamma_2

# Example usage with M = 10
M = 10
gamma_1, gamma_2 = find_optimal_steps(M)

print(f"For a condition number M = {M}:")
print(f"The optimal step sizes are:")
print(f"gamma_1 = {gamma_1}")
print(f"gamma_2 = {gamma_2}")

# Example usage with M = 100
M = 100
gamma_1, gamma_2 = find_optimal_steps(M)

print(f"\nFor a condition number M = {M}:")
print(f"The optimal step sizes are:")
print(f"gamma_1 = {gamma_1}")
print(f"gamma_2 = {gamma_2}")
