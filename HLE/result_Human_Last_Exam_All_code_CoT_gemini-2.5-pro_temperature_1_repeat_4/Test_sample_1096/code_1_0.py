import math

def solve_optimal_steps():
    """
    Calculates the optimal two step sizes for gradient descent for a given
    condition number M.
    """
    try:
        # The problem statement uses M to denote the condition number kappa.
        m_str = input("Enter the condition number M (kappa, must be >= 1): ")
        M = float(m_str)
        if M < 1:
            print("Error: Condition number M must be >= 1.")
            return
    except ValueError:
        print(f"Error: Invalid input '{m_str}'. Please enter a number.")
        return

    # Using the derived formulas for the optimal step sizes gamma_1 and gamma_2
    # Denominator D = M^2 + 6*M + 1
    denominator = M**2 + 6*M + 1

    # Common terms in the numerator
    term1 = 4 * (M + 1)
    term2 = 2 * math.sqrt(2) * (M - 1)

    # Calculate the two step sizes
    gamma1 = (term1 - term2) / denominator
    gamma2 = (term1 + term2) / denominator

    # Output the results, showing the numbers in the final equation
    print("\nBased on the Chebyshev polynomial solution for the two-step gradient method:")
    print(f"For a condition number M = {M}:")
    
    # Create formatted strings for the components of the formula
    s_term1 = f"4 * ({M} + 1)"
    s_term2 = f"2 * sqrt(2) * ({M} - 1)"
    s_denom = f"{M}**2 + 6 * {M} + 1"
    
    print("\nThe optimal step sizes are:")
    print(f"gamma_1 = ({s_term1} - {s_term2}) / ({s_denom})")
    print(f"gamma_1 = ({term1:.4f} - {term2:.4f}) / {denominator:.4f} = {gamma1}")
    
    print(f"\ngamma_2 = ({s_term1} + {s_term2}) / ({s_denom})")
    print(f"gamma_2 = ({term1:.4f} + {term2:.4f}) / {denominator:.4f} = {gamma2}")

if __name__ == '__main__':
    solve_optimal_steps()
