import math

def solve_optimal_steps():
    """
    This script calculates the optimal two-step gradient descent step sizes (gamma_1, gamma_2)
    for a given condition number M.
    """
    try:
        M_str = input("Enter the condition number M (kappa, M >= 1): ")
        M = float(M_str)
        if M < 1:
            print("M must be greater than or equal to 1.")
            return
    except ValueError:
        print("Invalid input. Please enter a number.")
        return

    # The optimal step sizes are derived from the roots of a polynomial related
    # to Chebyshev polynomials. The formulas for the two step sizes gamma_1 and gamma_2 are:
    # gamma_{1,2} = (4*(M+1) +/- 2*sqrt(2)*(M-1)) / (M^2 + 6*M + 1)
    
    # Let's print the formulas first
    print("\nThe optimal choice for the pair (gamma_1, gamma_2) is given by the formulas:")
    # Using unicode for better readability
    print("        4(M+1) \u2213 2\u221A2(M-1)")
    print(" \u03B3\u2081,\u2082 = \u2014\u2014\u2014\u2014\u2014\u2014\u2014\u2014\u2014\u2014\u2014")
    print("         M\u00B2 + 6M + 1")
    
    print("\nThese can be written out as two separate equations for each step size:")
    # Final equation parts
    num_m_plus_1 = 4
    num_m_minus_1_coeff = 2
    num_m_minus_1_sqrt = 2
    den_m_sq = 1
    den_m = 6
    den_const = 1
    
    # Equation 1
    print("\nEquation for gamma_1:")
    print(f"        {num_m_plus_1}(M+1) - {num_m_minus_1_coeff}\u221A{num_m_minus_1_sqrt}(M-1)")
    print(" \u03B3\u2081 = \u2014\u2014\u2014\u2014\u2014\u2014\u2014\u2014\u2014\u2014\u2014\u2014")
    print(f"        {den_m_sq if den_m_sq!=1 else ''}M\u00B2 + {den_m}M + {den_const}")
    
    # Equation 2
    print("\nEquation for gamma_2:")
    print(f"        {num_m_plus_1}(M+1) + {num_m_minus_1_coeff}\u221A{num_m_minus_1_sqrt}(M-1)")
    print(" \u03B3\u2082 = \u2014\u2014\u2014\u2014\u2014\u2014\u2014\u2014\u2014\u2014\u2014\u2014")
    print(f"        {den_m_sq if den_m_sq!=1 else ''}M\u00B2 + {den_m}M + {den_const}")
    
    # Calculate the numerical values for the given M
    denominator = M**2 + 6*M + 1
    sqrt_term = 2 * math.sqrt(2) * (M - 1)
    numerator_base = 4 * (M + 1)

    gamma1 = (numerator_base - sqrt_term) / denominator
    gamma2 = (numerator_base + sqrt_term) / denominator

    print(f"\nFor M = {M}:")
    print(f"gamma_1 = {gamma1}")
    print(f"gamma_2 = {gamma2}")

if __name__ == '__main__':
    solve_optimal_steps()