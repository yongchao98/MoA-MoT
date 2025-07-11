import math

def solve_optimal_steps():
    """
    Calculates the optimal two-step gradient descent learning rates (gammas)
    for a given condition number M.
    """
    # The condition number M (denoted as kappa in the prompt).
    # Since M is not specified in the problem, we will use a representative value.
    # A value of M=10 is used here as an example.
    # The user can change this value to their specific problem's M.
    M = 10.0

    print(f"Calculating optimal step sizes for a condition number M = {M}\n")

    # The formulas for the optimal gamma_1 and gamma_2 are derived from
    # a minimax problem involving Chebyshev polynomials. The resulting expressions are:
    # gamma_1,2 = (4*(M+1) +/- 2*sqrt(2)*(M-1)) / (M^2 + 6*M + 1)

    # Calculate the common denominator
    denominator = M**2 + 6*M + 1

    # Calculate the two terms in the numerator
    term1_num = 4 * (M + 1)
    term2_num = 2 * math.sqrt(2) * (M - 1)

    # Calculate the two step sizes. The order does not matter.
    gamma1 = (term1_num - term2_num) / denominator
    gamma2 = (term1_num + term2_num) / denominator
    
    # Let's ensure gamma1 is the smaller one for consistent output
    if gamma1 > gamma2:
        gamma1, gamma2 = gamma2, gamma1

    # Output the final equations with the numbers plugged in
    print("The optimal step sizes are the two roots of a quadratic equation derived from the Chebyshev polynomial.")
    print("The general formulas for the pair {gamma_1, gamma_2} are:")
    print("  (4(M+1) - 2*sqrt(2)*(M-1)) / (M^2+6M+1)   and   (4(M+1) + 2*sqrt(2)*(M-1)) / (M^2+6M+1)")
    print("\n--- Calculation Details ---\n")

    print("gamma_1 = (4 * ({M} + 1) - 2 * sqrt(2) * ({M} - 1)) / ({M}^2 + 6 * {M} + 1)".format(M=M))
    print("gamma_1 = (4 * {M_plus_1} - 2 * {sqrt2:.4f} * {M_minus_1}) / ({M_sq} + {M_6} + 1)".format(
        M_plus_1=M+1, sqrt2=math.sqrt(2), M_minus_1=M-1, M_sq=M**2, M_6=6*M))
    print("gamma_1 = ({term1:.4f} - {term2:.4f}) / {den:.4f}".format(term1=term1_num, term2=term2_num, den=denominator))
    print("gamma_1 = {num:.4f} / {den:.4f} = {result:.6f}".format(num=term1_num - term2_num, den=denominator, result=gamma1))
    
    print("\n")

    print("gamma_2 = (4 * ({M} + 1) + 2 * sqrt(2) * ({M} - 1)) / ({M}^2 + 6 * {M} + 1)".format(M=M))
    print("gamma_2 = (4 * {M_plus_1} + 2 * {sqrt2:.4f} * {M_minus_1}) / ({M_sq} + {M_6} + 1)".format(
        M_plus_1=M+1, sqrt2=math.sqrt(2), M_minus_1=M-1, M_sq=M**2, M_6=6*M))
    print("gamma_2 = ({term1:.4f} + {term2:.4f}) / {den:.4f}".format(term1=term1_num, term2=term2_num, den=denominator))
    print("gamma_2 = {num:.4f} / {den:.4f} = {result:.6f}".format(num=term1_num + term2_num, den=denominator, result=gamma2))

    print("\n--- Result ---")
    print("The best choice for the pair (gamma_1, gamma_2) is ({:.6f}, {:.6f}).".format(gamma1, gamma2))
    
    return gamma1, gamma2

if __name__ == '__main__':
    gamma1, gamma2 = solve_optimal_steps()
    # This format is for the final answer submission as requested by the prompt.
    # The result is for M=10.
    final_answer = f"<<<({gamma1:.6f}, {gamma2:.6f})>>>"
    # print(final_answer) # This print is commented out to not clutter the main output.

solve_optimal_steps()