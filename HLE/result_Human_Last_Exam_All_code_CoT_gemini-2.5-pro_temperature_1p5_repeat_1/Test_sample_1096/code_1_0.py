import math

def find_optimal_stepsizes():
    """
    Calculates the optimal two-step gradient descent learning rates (gamma_1, gamma_2)
    for a given condition number M.
    """
    # We will use a common example value for the condition number M.
    M = 10.0

    # The problem is to find gamma_1 and gamma_2 that minimize the worst-case
    # convergence rate ||x_2 - x_*||^2 / ||x_0 - x_*||^2.
    # This reduces to a minimax problem for a polynomial P(z) = (1-gamma_1*z)(1-gamma_2*z)
    # on the interval of scaled eigenvalues [1, M]. The solution is based on Chebyshev polynomials.

    if M <= 1:
        print("Condition number M must be greater than 1 for the problem to be non-trivial.")
        return

    # The derived optimal formulas for gamma_1 and gamma_2 are:
    # gamma_1 = (4 * (M + 1) + 2 * sqrt(2) * (M - 1)) / (M^2 + 6 * M + 1)
    # gamma_2 = (4 * (M + 1) - 2 * sqrt(2) * (M - 1)) / (M^2 + 6 * M + 1)

    print(f"### Calculating optimal step sizes for M = {M} ###")

    # Breaking down the formula to show each number
    M_plus_1 = M + 1
    M_minus_1 = M - 1
    denominator_val = M**2 + 6*M + 1
    sqrt_2 = math.sqrt(2)

    num_part1 = 4 * M_plus_1
    num_part2 = 2 * sqrt_2 * M_minus_1

    # Calculate gamma_1 and gamma_2
    gamma1_numerator = num_part1 + num_part2
    gamma1 = gamma1_numerator / denominator_val

    gamma2_numerator = num_part1 - num_part2
    gamma2 = gamma2_numerator / denominator_val
    
    print("\nThe general formulas for the optimal step sizes are:")
    print("gamma_1 = (4*(M + 1) + 2*sqrt(2)*(M - 1)) / (M^2 + 6*M + 1)")
    print("gamma_2 = (4*(M + 1) - 2*sqrt(2)*(M - 1)) / (M^2 + 6*M + 1)")

    print("\n--- Calculation with intermediate values ---")
    
    # Printing each number in the final equation
    print("\n[Denominator]:")
    print(f"M^2 + 6*M + 1 = {M:.1f}^2 + 6*{M:.1f} + 1 = {M**2} + {6*M} + 1 = {denominator_val}")

    print("\n[Numerator part 1]:")
    print(f"4 * (M + 1) = 4 * ({M:.1f} + 1) = 4 * {M_plus_1} = {num_part1}")
    
    print("\n[Numerator part 2]:")
    print(f"2 * sqrt(2) * (M - 1) = 2 * {sqrt_2:.4f} * ({M:.1f} - 1) = 2 * {sqrt_2:.4f} * {M_minus_1} = {num_part2:.4f}")

    print("\n--- Final Step-Size Values ---")

    print("\nFor gamma_1:")
    print(f"Numerator = {num_part1} + {num_part2:.4f} = {gamma1_numerator:.4f}")
    print(f"gamma_1 = Numerator / Denominator = {gamma1_numerator:.4f} / {denominator_val} = {gamma1:.6f}")

    print("\nFor gamma_2:")
    print(f"Numerator = {num_part1} - {num_part2:.4f} = {gamma2_numerator:.4f}")
    print(f"gamma_2 = Numerator / Denominator = {gamma2_numerator:.4f} / {denominator_val} = {gamma2:.6f}")

    print("\nNote: The quantity S=sqrt(M^2+(M-1)^2) provided in the problem description does not appear to be directly involved in the standard derivation of the optimal parameters.")
    
    # Per instructions, returning the final answer in the specified format
    # The "best choice" is the pair of these two values. The order does not affect the final convergence rate.
    print(f"\n<<<(gamma_1={gamma1:.6f}, gamma_2={gamma2:.6f})>>>")

if __name__ == '__main__':
    find_optimal_stepsizes()