import math

def calculate_optimal_steps(kappa):
    """
    Calculates the optimal two-step gradient descent learning rates.

    Args:
        kappa (float): The condition number of the function (M/m).

    Returns:
        tuple: A tuple containing the optimal (gamma_1, gamma_2).
    """
    # The optimal step sizes (gamma_1, gamma_2) are derived from the roots
    # of a shifted Chebyshev polynomial. The formulas are:
    # gamma_1 = (4*(kappa+1) - 2*sqrt(2)*(kappa-1)) / (kappa^2 + 6*kappa + 1)
    # gamma_2 = (4*(kappa+1) + 2*sqrt(2)*(kappa-1)) / (kappa^2 + 6*kappa + 1)

    # First, let's define the constants from the formulas.
    # The numbers in the final equation are: 4, 1, 2, 2, 1, 1, 2, 6, 1
    c_4 = 4.0
    c_2_sqrt_2 = 2.0 * math.sqrt(2)
    c_6 = 6.0
    c_1 = 1.0

    # Calculate the terms in the expressions
    kappa_plus_1 = kappa + c_1
    kappa_minus_1 = kappa - c_1

    # The denominator is the same for both gamma_1 and gamma_2
    denominator = kappa**2 + c_6 * kappa + c_1

    # Numerator for gamma_1
    numerator_1 = c_4 * kappa_plus_1 - c_2_sqrt_2 * kappa_minus_1
    # Numerator for gamma_2
    numerator_2 = c_4 * kappa_plus_1 + c_2_sqrt_2 * kappa_minus_1

    # Calculate the final values for gamma_1 and gamma_2
    gamma_1 = numerator_1 / denominator
    gamma_2 = numerator_2 / denominator

    # Print the results in a clear format
    print(f"For a condition number kappa = {kappa}:")
    print("\nThe optimal step sizes are:")
    print(f"gamma_1 = {gamma_1:.8f}")
    print(f"gamma_2 = {gamma_2:.8f}")
    print("\n--- Derivation details ---")
    print("Formula for gamma_1: (4 * (kappa + 1) - 2*sqrt(2) * (kappa - 1)) / (kappa^2 + 6*kappa + 1)")
    print(f"Calculation: ({c_4:.1f} * {kappa_plus_1:.1f} - {c_2_sqrt_2:.4f} * {kappa_minus_1:.1f}) / ({denominator:.1f}) = {gamma_1:.8f}")
    print("\nFormula for gamma_2: (4 * (kappa + 1) + 2*sqrt(2) * (kappa - 1)) / (kappa^2 + 6*kappa + 1)")
    print(f"Calculation: ({c_4:.1f} * {kappa_plus_1:.1f} + {c_2_sqrt_2:.4f} * {kappa_minus_1:.1f}) / ({denominator:.1f}) = {gamma_2:.8f}")

if __name__ == '__main__':
    # You can change this value to compute the step sizes for a different condition number.
    condition_number_kappa = 10.0
    calculate_optimal_steps(condition_number_kappa)