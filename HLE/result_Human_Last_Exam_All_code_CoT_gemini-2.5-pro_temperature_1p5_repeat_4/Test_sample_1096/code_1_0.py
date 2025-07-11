import math

def get_optimal_gammas(M):
    """
    Computes the optimal step sizes (gamma_1, gamma_2) for a 2-step gradient descent.

    Args:
        M (float): The condition number kappa. Assumed to be > 1.

    Returns:
        tuple: A pair of floats (gamma_1, gamma_2).
    """
    if M <= 1:
        raise ValueError("Condition number M must be greater than 1.")

    # The formulas are derived from Chebyshev polynomial properties.
    # gamma_1,2 = (4*(M+1) +- 2*sqrt(2)*(M-1)) / (M^2 + 6*M + 1)
    
    kappa = M
    
    print("The formulas for the optimal pair (gamma_1, gamma_2) are:")
    print("gamma_1 = (4*(kappa+1) - 2*sqrt(2)*(kappa-1)) / (kappa^2 + 6*kappa + 1)")
    print("gamma_2 = (4*(kappa+1) + 2*sqrt(2)*(kappa-1)) / (kappa^2 + 6*kappa + 1)\n")
    
    # As requested, printing the numbers in the final equations for the formulas above.
    # The constants are 4, 1, 2, 2, 1 (for numerator) and 2, 6, 1 (for denominator).
    # Since sqrt(2) involves the number 2, it is included.
    print("The constant numbers appearing in the final equations are:")
    print("Numerator constants: 4, 1, 2")
    print("Denominator constants: 6, 1\n")


    denominator = kappa**2 + 6*kappa + 1
    term1 = 4 * (kappa + 1)
    term2 = 2 * math.sqrt(2) * (kappa - 1)

    gamma1 = (term1 - term2) / denominator
    gamma2 = (term1 + term2) / denominator

    return gamma1, gamma2

def main():
    """
    Main function to demonstrate the calculation for a sample value of M.
    """
    # Let's use a sample value for the condition number, e.g., M = 10
    M = 10.0
    
    print(f"Calculating optimal gammas for M = {M}:\n")
    
    gamma1, gamma2 = get_optimal_gammas(M)
    
    print(f"Numerical values for M = {M}:")
    print(f"gamma_1 = {gamma1}")
    print(f"gamma_2 = {gamma2}")
    
    # The final answer in the requested format
    # The question asks for the best choice of the pair, which are the formulas.
    answer_string = "gamma_1,2 = (4*(M+1) +- 2*sqrt(2)*(M-1)) / (M^2 + 6*M + 1)"
    # print(f"\n<<<{answer_string}>>>")


if __name__ == "__main__":
    main()
