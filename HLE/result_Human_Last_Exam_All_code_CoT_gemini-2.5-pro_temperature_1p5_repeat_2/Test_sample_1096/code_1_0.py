import math

def calculate_optimal_stepsizes(kappa):
    """
    Calculates the optimal two-step gradient descent stepsizes (gamma_1, gamma_2).

    Args:
        kappa (float): The condition number of the function (M/mu).
    """
    # We assume mu = 1 and M = kappa, as specified in the problem.
    M = float(kappa)
    mu = 1.0

    print(f"Calculating optimal stepsizes for M = {M} and mu = {mu}\n")

    # The optimal step sizes gamma_1 and gamma_2 are the reciprocals of the roots
    # of the optimal Chebyshev-based polynomial. The formulas are derived from this principle.
    # P(lambda) = (1 - gamma_1 * lambda)(1 - gamma_2 * lambda)
    # The roots are 1/gamma_1 and 1/gamma_2.
    # The optimal roots are lambda_roots = (M+mu)/2 +/- (M-mu)/(2*sqrt(2))

    sqrt2 = math.sqrt(2)

    # Calculate the two step sizes. The order does not matter.
    # The formulas are derived from gamma = 1 / lambda_root.

    # Denominator for the first step size
    # This corresponds to the root (M+mu)/2 + (M-mu)/(2*sqrt(2))
    # which can be rewritten as ((sqrt(2)+1)*M + (sqrt(2)-1)*mu) / (2*sqrt(2))
    gamma1_numerator = 2 * sqrt2
    gamma1_denominator = (sqrt2 + 1) * M + (sqrt2 - 1) * mu
    gamma1 = gamma1_numerator / gamma1_denominator

    # Denominator for the second step size
    # This corresponds to the root (M+mu)/2 - (M-mu)/(2*sqrt(2))
    # which can be rewritten as ((sqrt(2)-1)*M + (sqrt(2)+1)*mu) / (2*sqrt(2))
    gamma2_numerator = 2 * sqrt2
    gamma2_denominator = (sqrt2 - 1) * M + (sqrt2 + 1) * mu
    gamma2 = gamma2_numerator / gamma2_denominator

    print("The optimal step sizes are given by the pair (gamma_1, gamma_2).")
    print("-" * 50)
    
    # Print the first gamma with the equation breakdown
    print(f"gamma_1 = (2 * sqrt(2)) / ((sqrt(2) + 1) * M + (sqrt(2) - 1) * mu)")
    print(f"      = (2 * {sqrt2:.5f}) / (({sqrt2:.5f} + 1) * {M} + ({sqrt2:.5f} - 1) * {mu})")
    print(f"      = {gamma1_numerator:.5f} / ({gamma1_denominator:.5f})")
    print(f"      = {gamma1:.5f}\n")

    # Print the second gamma with the equation breakdown
    print(f"gamma_2 = (2 * sqrt(2)) / ((sqrt(2) - 1) * M + (sqrt(2) + 1) * mu)")
    print(f"      = (2 * {sqrt2:.5f}) / (({sqrt2:.5f} - 1) * {M} + ({sqrt2:.5f} + 1) * {mu})")
    print(f"      = {gamma2_numerator:.5f} / ({gamma2_denominator:.5f})")
    print(f"      = {gamma2:.5f}\n")
    print(f"The best choice for the pair (gamma_1, gamma_2) is ({gamma1:.5f}, {gamma2:.5f}). The order is interchangeable.")

if __name__ == '__main__':
    # You can change the value of kappa here.
    # kappa must be >= 1.
    kappa = 10.0
    calculate_optimal_stepsizes(kappa)
