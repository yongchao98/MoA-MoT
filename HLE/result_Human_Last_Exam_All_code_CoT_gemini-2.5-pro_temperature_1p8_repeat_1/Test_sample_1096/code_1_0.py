import math

def calculate_optimal_steps():
    """
    Calculates the optimal two-step gradient descent learning rates
    for a given condition number M.
    """
    # We use M as the condition number, as specified in the problem (M=kappa)
    # Let's use a sample value for M for demonstration purposes.
    M = 10.0

    print(f"This script calculates the optimal step sizes (gamma_1, gamma_2) for M = {M}\n")

    # The formulas for the optimal step sizes are derived from Chebyshev polynomials:
    # gamma_1,2 = (4*(M+1) +/- 2*sqrt(2)*(M-1)) / (M^2 + 6*M + 1)
    # The order of gamma_1 and gamma_2 (which one has + or -) is arbitrary.
    # We will assign the smaller step to gamma_1.

    # Calculate the individual terms for clarity
    numerator_term1 = 4 * (M + 1)
    numerator_term2 = 2 * math.sqrt(2) * (M - 1)
    denominator = M**2 + 6*M + 1

    # Calculate gamma_1 and gamma_2
    gamma1 = (numerator_term1 - numerator_term2) / denominator
    gamma2 = (numerator_term1 + numerator_term2) / denominator

    # The user asked to output each number in the final equation.
    print("The general formulas for the optimal step sizes are:")
    print("gamma_1 = (4*(M+1) - 2*sqrt(2)*(M-1)) / (M^2 + 6*M + 1)")
    print("gamma_2 = (4*(M+1) + 2*sqrt(2)*(M-1)) / (M^2 + 6*M + 1)\n")
    
    # --- Output for gamma_1 ---
    print(f"Calculation for gamma_1 (smaller step) with M = {M}:")
    # Step 1: Show the formula with M substituted
    print(f"gamma_1 = (4*({M} + 1) - 2*sqrt(2)*({M} - 1)) / ({M}^2 + 6*{M} + 1)")
    # Step 2: Show the values of the components
    print(f"gamma_1 = ({numerator_term1} - 2 * {math.sqrt(2):.5f} * {M - 1}) / ({M**2} + {6*M} + 1)")
    # Step 3: Show the calculated terms
    print(f"gamma_1 = ({numerator_term1} - {numerator_term2:.5f}) / {denominator}")
    # Step 4: Show the final result
    print(f"gamma_1 = {gamma1:.6f}\n")
    
    # --- Output for gamma_2 ---
    print(f"Calculation for gamma_2 (larger step) with M = {M}:")
    # Step 1: Show the formula with M substituted
    print(f"gamma_2 = (4*({M} + 1) + 2*sqrt(2)*({M} - 1)) / ({M}^2 + 6*{M} + 1)")
    # Step 2: Show the values of the components
    print(f"gamma_2 = ({numerator_term1} + 2 * {math.sqrt(2):.5f} * {M - 1}) / ({M**2} + {6*M} + 1)")
    # Step 3: Show the calculated terms
    print(f"gamma_2 = ({numerator_term1} + {numerator_term2:.5f}) / {denominator}")
    # Step 4: Show the final result
    print(f"gamma_2 = {gamma2:.6f}")


if __name__ == '__main__':
    calculate_optimal_steps()