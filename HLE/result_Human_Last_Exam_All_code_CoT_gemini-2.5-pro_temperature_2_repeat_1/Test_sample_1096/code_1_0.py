import math

def calculate_optimal_steps():
    """
    Calculates the optimal two-step gradient descent learning rates (gammas)
    for a given condition number M.
    """
    try:
        m_str = input("Enter the condition number M (must be > 1): ")
        M = float(m_str)
        if M <= 1:
            print("Error: The condition number M must be greater than 1.")
            return
    except ValueError:
        print(f"Error: Invalid input '{m_str}'. Please enter a number.")
        return

    # Based on the Chebyshev polynomial solution, the optimal step sizes gamma_1 and gamma_2
    # are given by the following formulas:
    # gamma = (4*(M+1) +/- 2*sqrt(2)*(M-1)) / (M^2 + 6*M + 1)
    
    # Let's compute each component of the final equations.
    # Let the formulas be:
    # gamma_1 = (N1 - N2) / D
    # gamma_2 = (N1 + N2) / D
    #
    # where:
    # N1 = 4 * (M+1)
    # N2 = 2 * sqrt(2) * (M-1)
    # D = M^2 + 6*M + 1

    n1 = 4 * (M + 1)
    n2 = 2 * math.sqrt(2) * (M - 1)
    d = M**2 + 6 * M + 1

    # Now compute gamma_1 and gamma_2
    # By convention, gamma_1 is the smaller step and gamma_2 is the larger one.
    gamma_1 = (n1 - n2) / d
    gamma_2 = (n1 + n2) / d

    # Output the equations with the calculated numbers
    print("\n" + "="*50)
    print(f"Optimal step sizes for M = {M}:")
    print("="*50)
    print("The general formulas for the optimal step sizes are:")
    print("  gamma_1 = (4*(M+1) - 2*sqrt(2)*(M-1)) / (M^2 + 6*M + 1)")
    print("  gamma_2 = (4*(M+1) + 2*sqrt(2)*(M-1)) / (M^2 + 6*M + 1)")
    
    print("\n" + "-"*50)
    print(f"Substituting M = {M} into the formulas:")
    print("-" * 50)
    
    # As requested, printing each number in the final equation.
    print("\nCalculation for gamma_1 (the smaller step size):")
    print(f"  Numerator term 1         (4*(M+1)): {n1:.4f}")
    print(f"  Numerator term 2  (2*sqrt(2)*(M-1)): {n2:.4f}")
    print(f"  Denominator        (M^2 + 6*M + 1): {d:.4f}")
    print("-" * 50)
    print(f"  gamma_1 = ({n1:.4f} - {n2:.4f}) / {d:.4f}")
    print(f"  Final value for gamma_1 = {gamma_1:.6f}")

    print("\nCalculation for gamma_2 (the larger step size):")
    print(f"  Numerator term 1         (4*(M+1)): {n1:.4f}")
    print(f"  Numerator term 2  (2*sqrt(2)*(M-1)): {n2:.4f}")
    print(f"  Denominator        (M^2 + 6*M + 1): {d:.4f}")
    print("-" * 50)
    print(f"  gamma_2 = ({n1:.4f} + {n2:.4f}) / {d:.4f}")
    print(f"  Final value for gamma_2 = {gamma_2:.6f}")
    print("="*50)

if __name__ == '__main__':
    calculate_optimal_steps()