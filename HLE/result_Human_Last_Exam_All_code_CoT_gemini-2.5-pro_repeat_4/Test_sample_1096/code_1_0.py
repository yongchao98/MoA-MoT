import math

def calculate_optimal_steps(M):
    """
    Calculates the optimal two-step gradient descent learning rates.

    Args:
        M (float): The condition number kappa of the function.

    Returns:
        tuple: A tuple containing the two optimal learning rates (gamma_1, gamma_2).
    """
    if M <= 1:
        print("Condition number M must be greater than 1.")
        return None, None
        
    print(f"Calculating for M = {M}")

    # Denominator of the expressions for the learning rates
    denominator = M**2 + 6*M + 1
    print(f"Denominator (M^2 + 6M + 1) = {denominator}")

    # Numerator terms
    term1_val = 4 * (M + 1)
    term2_val = 2 * math.sqrt(2) * (M - 1)
    print(f"Numerator term 1 (4 * (M + 1)) = {term1_val}")
    print(f"Numerator term 2 (2 * sqrt(2) * (M - 1)) = {term2_val:.4f}")
    
    # Calculate the two learning rates
    gamma_1 = (term1_val - term2_val) / denominator
    gamma_2 = (term1_val + term2_val) / denominator

    return gamma_1, gamma_2

# --- Main execution ---
# Set the condition number M (kappa)
# For example, let's use a common value M = 10
M = 10

# Calculate and print the optimal step sizes
gamma_1, gamma_2 = calculate_optimal_steps(M)

if gamma_1 is not None:
    print("\n--- Optimal Step Sizes ---")
    print(f"The equation for gamma_1 is: (4*({M}+1) - 2*sqrt(2)*({M}-1)) / ({M}^2+6*{M}+1)")
    print(f"gamma_1 = {gamma_1:.4f}")
    print(f"The equation for gamma_2 is: (4*({M}+1) + 2*sqrt(2)*({M}-1)) / ({M}^2+6*{M}+1)")
    print(f"gamma_2 = {gamma_2:.4f}")
