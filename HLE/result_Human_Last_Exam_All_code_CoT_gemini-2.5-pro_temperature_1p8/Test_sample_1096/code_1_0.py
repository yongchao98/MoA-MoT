import math

def calculate_optimal_stepsizes(M):
    """
    Calculates the optimal step sizes for a two-step gradient descent method.

    Args:
        M (float): The condition number kappa. Must be >= 1.
    """
    if M < 1:
        print("Error: M (kappa) must be greater than or equal to 1.")
        return

    # Numerator parts
    # The term 4*(M+1)
    num_part1 = 4 * (M + 1)
    # The term 2*sqrt(2)*(M-1)
    num_part2 = 2 * math.sqrt(2) * (M - 1)
    
    # Denominator
    # The term M^2 + 6*M + 1
    denominator = M**2 + 6 * M + 1

    # Calculate the two step sizes
    gamma1 = (num_part1 - num_part2) / denominator
    gamma2 = (num_part1 + num_part2) / denominator
    
    print("For a condition number M = {}".format(M))
    print("\n--- The formula for the optimal step sizes gamma_1,2 is: (A +/- B) / C ---")
    print("A = 4 * (M + 1) = {:.4f}".format(num_part1))
    print("B = 2 * sqrt(2) * (M - 1) = {:.4f}".format(num_part2))
    print("C = M^2 + 6*M + 1 = {:.4f}".format(denominator))
    print("\nThe optimal step sizes are:")
    print("gamma_1 = {:.8f}".format(gamma1))
    print("gamma_2 = {:.8f}".format(gamma2))

# Example usage with a condition number M=10
# You can change this value to test with other condition numbers.
M = 10
calculate_optimal_stepsizes(M)

# Final answer format for the derived formulas.
# Note: The order of gamma_1 and gamma_2 can be swapped.
# <<<gamma_1 = (4*(M+1) - 2*sqrt(2)*(M-1))/(M**2 + 6*M + 1), gamma_2 = (4*(M+1) + 2*sqrt(2)*(M-1))/(M**2 + 6*M + 1)>>>