import math

def get_optimal_step_sizes(M):
    """
    Calculates the optimal two-step gradient descent step sizes (gamma_1, gamma_2)
    for an M-smooth and 1-strongly convex function.

    Args:
        M (float): The condition number (smoothness / strong convexity).

    Returns:
        tuple: A tuple containing the two optimal step sizes (gamma_1, gamma_2).
    """
    if M < 1:
        raise ValueError("Condition number M must be >= 1")
    if M == 1:
        # For M=1, the function is a simple quadratic with a single eigenvalue.
        # The optimal step is 1/M = 1. The formula gives 8/8=1.
        return (1.0, 1.0)

    # Denominator of the expression for gamma_1 and gamma_2
    denominator = M**2 + 6*M + 1
    
    # Numerator parts
    term1 = 4 * (M + 1)
    term2 = 2 * math.sqrt(2) * (M - 1)
    
    # Calculate gamma_1 and gamma_2
    gamma_1 = (term1 - term2) / denominator
    gamma_2 = (term1 + term2) / denominator
    
    return gamma_1, gamma_2

def main():
    """
    Main function to execute the calculation and print the results.
    """
    # The condition number M=kappa. You can change this value.
    # M=1 is a trivial case, M > 1 is interesting.
    M = 10.0
    
    print(f"For a function with condition number M = {M}:")
    
    gamma_1, gamma_2 = get_optimal_step_sizes(M)
    
    # Print the equation with numbers
    # We show the final values in the equations for clarity.
    
    # Numerator of gamma_1
    num1 = 4 * (M + 1) - 2 * math.sqrt(2) * (M - 1)
    # Numerator of gamma_2
    num2 = 4 * (M + 1) + 2 * math.sqrt(2) * (M - 1)
    # Denominator is common
    den = M**2 + 6*M + 1

    print("\nThe optimal step sizes (gamma_1, gamma_2) are calculated as:")
    print(f"gamma_1 = (4*({M}+1) - 2*sqrt(2)*({M}-1)) / ({M}^2 + 6*{M} + 1)")
    print(f"        = ({4*(M+1):.4f} - {2*math.sqrt(2)*(M-1):.4f}) / {den:.4f}")
    print(f"        = {num1:.4f} / {den:.4f}")
    print(f"        = {gamma_1:.10f}")
    
    print(f"\ngamma_2 = (4*({M}+1) + 2*sqrt(2)*({M}-1)) / ({M}^2 + 6*{M} + 1)")
    print(f"        = ({4*(M+1):.4f} + {2*math.sqrt(2)*(M-1):.4f}) / {den:.4f}")
    print(f"        = {num2:.4f} / {den:.4f}")
    print(f"        = {gamma_2:.10f}")

    print("\n----------------------------------------------------")
    print("Final calculated pair (gamma_1, gamma_2):")
    print((gamma_1, gamma_2))
    print("----------------------------------------------------")

if __name__ == "__main__":
    main()