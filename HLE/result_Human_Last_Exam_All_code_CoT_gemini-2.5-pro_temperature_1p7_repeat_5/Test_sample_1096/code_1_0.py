import math

def compute_optimal_steps(M):
    """
    Computes the optimal pair of step sizes (gamma_1, gamma_2) for a two-step
    gradient descent method.

    Args:
        M (float): The smoothness constant (conditioning number kappa, assuming mu=1).
                   It must be greater than 1.
    """
    if M <= 1:
        print("M (kappa) must be greater than 1 for the problem to be non-trivial.")
        return

    kappa = M

    # The optimal step sizes are derived from the coefficients of the optimal
    # polynomial P(lambda) that minimizes the max norm on [1, kappa].
    # P(lambda) = 1 - (gamma_1 + gamma_2) * lambda + (gamma_1 * gamma_2) * lambda^2
    # The coefficients are determined by the Chebyshev polynomial of degree 2.
    
    # Denominator for the expressions of gamma_1 and gamma_2
    denominator = kappa**2 + 6*kappa + 1
    
    # The sum and product of the optimal step sizes
    gamma_sum = 8 * (kappa + 1) / denominator
    gamma_prod = 8 / denominator
    
    # We solve for gamma_1 and gamma_2 from their sum and product.
    # They are the roots of the quadratic equation: x^2 - gamma_sum * x + gamma_prod = 0
    discriminant_sqrt = math.sqrt(gamma_sum**2 - 4 * gamma_prod)
    
    gamma_1 = (gamma_sum - discriminant_sqrt) / 2
    gamma_2 = (gamma_sum + discriminant_sqrt) / 2
    
    # Alternative direct formulas:
    # num_part1 = 4 * (kappa + 1)
    # num_part2 = 2 * math.sqrt(2) * (kappa - 1)
    # gamma_1_alt = (num_part1 - num_part2) / denominator
    # gamma_2_alt = (num_part1 + num_part2) / denominator
    
    print(f"For M = {M}:")
    print("----------------")
    print(f"Denominator value (kappa^2 + 6*kappa + 1) = {denominator}")
    print(f"Numerator term 1 (4 * (kappa + 1)) = {4 * (kappa + 1)}")
    print(f"Numerator term 2 (2 * sqrt(2) * (kappa - 1)) = {2 * math.sqrt(2) * (kappa - 1)}")
    
    # Printing the final equations with calculated values
    print("\nThe optimal step sizes (gamma_1, gamma_2) are:")
    print(f"gamma_1 = ({4 * (kappa + 1)} - {2 * math.sqrt(2) * (kappa - 1)}) / {denominator} = {gamma_1}")
    print(f"gamma_2 = ({4 * (kappa + 1)} + {2 * math.sqrt(2) * (kappa - 1)}) / {denominator} = {gamma_2}")


if __name__ == '__main__':
    # You can change the value of M (the smoothness constant) here.
    # It must be a value > 1. A common example value is 10.
    M_val = 10.0
    compute_optimal_steps(M_val)
