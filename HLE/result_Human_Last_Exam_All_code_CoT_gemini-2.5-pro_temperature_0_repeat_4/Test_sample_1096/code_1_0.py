import math

def find_optimal_stepsizes():
    """
    Calculates the optimal two-step gradient descent stepsizes (gamma1, gamma2)
    for a given condition number M.
    """
    # Set the condition number M. You can change this value.
    M = 10.0

    # The problem is to find the pair (gamma1, gamma2) that minimizes
    # the convergence rate. The solution is derived using Chebyshev polynomials.

    # The formulas for the optimal step sizes gamma1 and gamma2 are:
    # gamma1 = (4*(M+1) - 2*sqrt(2)*(M-1)) / (M^2+6*M+1)
    # gamma2 = (4*(M+1) + 2*sqrt(2)*(M-1)) / (M^2+6*M+1)

    print(f"Calculating optimal step sizes for M = {M}\n")

    # Breaking down the formula into parts
    M_plus_1 = M + 1
    M_minus_1 = M - 1
    denominator = M**2 + 6*M + 1
    sqrt_2 = math.sqrt(2)

    # Calculate the numerators for gamma1 and gamma2
    gamma1_numerator = 4 * M_plus_1 - 2 * sqrt_2 * M_minus_1
    gamma2_numerator = 4 * M_plus_1 + 2 * sqrt_2 * M_minus_1

    # Calculate the final values for gamma1 and gamma2
    gamma1 = gamma1_numerator / denominator
    gamma2 = gamma2_numerator / denominator

    # Output the equations with the numbers plugged in
    print("Formula for gamma1:")
    print(f"gamma1 = (4 * ({M} + 1) - 2 * sqrt(2) * ({M} - 1)) / ({M}^2 + 6 * {M} + 1)")
    print(f"gamma1 = (4 * {M_plus_1} - 2 * {sqrt_2:.4f} * {M_minus_1}) / {denominator}")
    print(f"gamma1 = ({4 * M_plus_1} - {2 * sqrt_2 * M_minus_1:.4f}) / {denominator}")
    print(f"gamma1 = {gamma1_numerator:.4f} / {denominator}")
    print(f"gamma1 = {gamma1:.6f}\n")

    print("Formula for gamma2:")
    print(f"gamma2 = (4 * ({M} + 1) + 2 * sqrt(2) * ({M} - 1)) / ({M}^2 + 6 * {M} + 1)")
    print(f"gamma2 = (4 * {M_plus_1} + 2 * {sqrt_2:.4f} * {M_minus_1}) / {denominator}")
    print(f"gamma2 = ({4 * M_plus_1} + {2 * sqrt_2 * M_minus_1:.4f}) / {denominator}")
    print(f"gamma2 = {gamma2_numerator:.4f} / {denominator}")
    print(f"gamma2 = {gamma2:.6f}\n")

    print(f"The best choice for the pair (gamma1, gamma2) is ({gamma1:.6f}, {gamma2:.6f}).")

if __name__ == '__main__':
    find_optimal_stepsizes()