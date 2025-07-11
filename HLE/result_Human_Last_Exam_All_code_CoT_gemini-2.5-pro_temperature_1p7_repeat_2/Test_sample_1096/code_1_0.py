import math

def calculate_optimal_steps(kappa):
    """
    Calculates the optimal two-step gradient descent learning rates (gammas).

    The function f is M-smooth and mu-strongly convex.
    The algorithm is:
    x_1 = x_0 - gamma_1 * grad(f(x_0))
    x_2 = x_1 - gamma_2 * grad(f(x_1))

    We assume mu=1 and M=kappa (the condition number).
    This function finds the pair (gamma_1, gamma_2) that minimizes
    the uniform convergence rate ||x_2 - x*||^2 / ||x_0 - x*||^2.

    Args:
        kappa (float): The condition number M/mu, where mu=1. Must be > 1.
    """
    if not isinstance(kappa, (int, float)) or kappa <= 1:
        print("Error: Kappa must be a number greater than 1.")
        return

    # The optimal step sizes are the roots of the quadratic equation x^2 - s*x + p = 0
    # where s and p are derived from the Chebyshev polynomial approximation.
    # The formulas are derived as:
    # gamma_1,2 = (4*(kappa+1) +/- 2*sqrt(2)*(kappa-1)) / (kappa^2+6*kappa+1)

    # Denominator
    denominator = kappa**2 + 6*kappa + 1

    # Numerator terms
    num_term1 = 4 * (kappa + 1)
    num_term2 = 2 * math.sqrt(2) * (kappa - 1)

    # Calculate gamma_1 and gamma_2
    gamma1 = (num_term1 - num_term2) / denominator
    gamma2 = (num_term1 + num_term2) / denominator

    # Print the results in a structured way
    print(f"For a condition number kappa = {kappa}:")
    print("-" * 40)
    print("The optimal step sizes (gamma_1, gamma_2) are calculated using the formula:")
    print("(4*(kappa+1) \u00B1 2*\u221A2*(kappa-1)) / (kappa\u00B2 + 6*kappa + 1)\n")
    
    print("Let's break down the calculation:")
    
    # Printing each number in the final equation
    print("Numerator part 1 (4*(kappa+1)):")
    print(f"4 * ({kappa} + 1) = {num_term1}")

    print("\nNumerator part 2 (2*\u221A2*(kappa-1)):")
    print(f"2 * {math.sqrt(2):.4f} * ({kappa} - 1) = {num_term2:.4f}")
    
    print("\nDenominator (kappa\u00B2 + 6*kappa + 1):")
    print(f"{kappa}\u00B2 + 6*{kappa} + 1 = {denominator}\n")

    print("-" * 40)
    print("Final optimal step sizes:")
    print(f"\u03B3\u2081 = ({num_term1:.4f} - {num_term2:.4f}) / {denominator:.4f} = {gamma1:.6f}")
    print(f"\u03B3\u2082 = ({num_term1:.4f} + {num_term2:.4f}) / {denominator:.4f} = {gamma2:.6f}")
    print("-" * 40)

if __name__ == '__main__':
    # You can change this value to compute the steps for a different kappa
    try:
        kappa_input = float(input("Enter the condition number kappa (e.g., 10, 100): "))
        calculate_optimal_steps(kappa_input)
    except ValueError:
        print("Invalid input. Please enter a number.")
    except Exception as e:
        print(f"An error occurred: {e}")
