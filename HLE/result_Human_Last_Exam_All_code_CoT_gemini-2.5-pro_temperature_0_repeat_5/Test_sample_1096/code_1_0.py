import math

def calculate_optimal_step_sizes(M):
    """
    Calculates the optimal two-step gradient descent step sizes (gamma_1, gamma_2)
    for M-smooth and 1-strongly convex functions.

    Args:
        M (float): The smoothness constant (and condition number since m=1).

    Returns:
        tuple: A tuple containing the optimal step sizes (gamma_1, gamma_2).
    """
    if M <= 1:
        print("M (kappa) must be greater than 1.")
        return None, None

    # The problem is to find gamma_1, gamma_2 that minimize max_{lambda in [1, M]} |(1-gamma_1*lambda)(1-gamma_2*lambda)|
    # The solution is based on Chebyshev polynomials. The optimal step sizes are the
    # reciprocals of the roots of the scaled Chebyshev polynomial of degree 2.
    # The formulas for the step sizes are:
    # gamma_1 = (4(M+1) - 2*sqrt(2)*(M-1)) / (M^2 + 6*M + 1)
    # gamma_2 = (4(M+1) + 2*sqrt(2)*(M-1)) / (M^2 + 6*M + 1)

    # Let's calculate the numerator and denominator for clarity.
    # Denominator is common for both
    denominator = M**2 + 6*M + 1

    # Numerator for gamma_1
    numerator1 = 4 * (M + 1) - 2 * math.sqrt(2) * (M - 1)

    # Numerator for gamma_2
    numerator2 = 4 * (M + 1) + 2 * math.sqrt(2) * (M - 1)

    # Calculate gamma_1 and gamma_2
    gamma1 = numerator1 / denominator
    gamma2 = numerator2 / denominator

    print(f"For M = {M}:")
    print("\n--- Calculation for gamma_1 ---")
    print(f"The equation is: (4*(M+1) - 2*sqrt(2)*(M-1)) / (M^2 + 6*M + 1)")
    print(f"Numerator: 4*({M}+1) - 2*sqrt(2)*({M}-1) = {numerator1}")
    print(f"Denominator: {M}^2 + 6*{M} + 1 = {denominator}")
    print(f"gamma_1 = {numerator1} / {denominator} = {gamma1}")

    print("\n--- Calculation for gamma_2 ---")
    print(f"The equation is: (4*(M+1) + 2*sqrt(2)*(M-1)) / (M^2 + 6*M + 1)")
    print(f"Numerator: 4*({M}+1) + 2*sqrt(2)*({M}-1) = {numerator2}")
    print(f"Denominator: {M}^2 + 6*{M} + 1 = {denominator}")
    print(f"gamma_2 = {numerator2} / {denominator} = {gamma2}")

    return gamma1, gamma2

if __name__ == '__main__':
    # Example usage with M = 10 (a common example for condition number)
    # The user can change this value.
    M_kappa = 10.0
    calculate_optimal_step_sizes(M_kappa)
