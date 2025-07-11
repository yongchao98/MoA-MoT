import math

def calculate_optimal_steps(M):
    """
    Calculates the optimal two-step gradient descent step sizes (gamma_1, gamma_2)
    for an M-smooth and 1-strongly convex function.

    Args:
        M (float): The conditioning number (kappa), M > 1.
    """
    if M <= 1:
        print("M (kappa) must be greater than 1.")
        return

    print(f"Calculating optimal step sizes for M = {M}:\n")

    # The formulas for gamma_1 and gamma_2 are derived from the coefficients
    # of the optimal Chebyshev polynomial for this problem.
    # gamma_1, gamma_2 are the roots of the quadratic equation:
    # (M^2 + 6M + 1) * z^2 - 8(M+1) * z + 8 = 0
    # The solution is z = (4(M+1) +/- 2*sqrt(2)*(M-1)) / (M^2 + 6M + 1)

    # Let's calculate the components of this formula.
    # This corresponds to the request to "output each number in the final equation".
    numerator_term1 = 4 * (M + 1)
    numerator_term2 = 2 * math.sqrt(2) * (M - 1)
    denominator = M**2 + 6 * M + 1

    print("The final equation for the step sizes gamma is:")
    print("gamma = (numerator_term1 +/- numerator_term2) / denominator\n")
    print(f"where:")
    print(f"  numerator_term1 = 4 * (M + 1) = {numerator_term1:.4f}")
    print(f"  numerator_term2 = 2 * sqrt(2) * (M - 1) = {numerator_term2:.4f}")
    print(f"  denominator = M^2 + 6*M + 1 = {denominator:.4f}\n")


    # Calculate the two optimal step sizes
    gamma_1 = (numerator_term1 + numerator_term2) / denominator
    gamma_2 = (numerator_term1 - numerator_term2) / denominator

    print("The optimal step sizes (gamma_1, gamma_2) are:")
    print(f"gamma_1 = ({numerator_term1:.4f} + {numerator_term2:.4f}) / {denominator:.4f} = {gamma_1:.4f}")
    print(f"gamma_2 = ({numerator_term1:.4f} - {numerator_term2:.4f}) / {denominator:.4f} = {gamma_2:.4f}")
    print("\nNote: The assignment to gamma_1 and gamma_2 is interchangeable.")


# --- User Execution ---
# You can change the value of M (the conditioning number) here.
# For example, let's use M = 10.
M = 10
calculate_optimal_steps(M)