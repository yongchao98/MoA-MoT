import math

def find_optimal_steps(kappa):
    """
    Calculates the optimal two-step gradient descent learning rates (gamma_1, gamma_2).

    The problem is to minimize the worst-case convergence rate for M-smooth and
    m-strongly convex functions, where kappa = M/m is the condition number.
    We assume m=1 and M=kappa.

    The optimal step sizes are the roots of the quadratic equation:
    (kappa^2 + 6*kappa + 1) * z^2 - 8*(kappa+1)*z + 8 = 0

    The explicit formulas for the roots gamma_1 and gamma_2 are:
    gamma_1 = (4*(kappa+1) - 2*sqrt(2)*(kappa-1)) / (kappa^2 + 6*kappa + 1)
    gamma_2 = (4*(kappa+1) + 2*sqrt(2)*(kappa-1)) / (kappa^2 + 6*kappa + 1)

    Args:
        kappa (float): The condition number of the function (M/m). Must be >= 1.

    Returns:
        tuple: A tuple containing the optimal step sizes (gamma_1, gamma_2).
    """
    if kappa < 1:
        raise ValueError("Condition number kappa must be >= 1.")
    if kappa == 1:
        # For kappa=1, the function is perfectly conditioned.
        # The formula gives gamma1=gamma2=1, which is the exact step to the minimum.
        return (1.0, 1.0)

    # Denominator of the expressions for gamma_1 and gamma_2
    denominator = kappa**2 + 6 * kappa + 1

    # Numerator parts
    term1 = 4 * (kappa + 1)
    term2 = 2 * math.sqrt(2) * (kappa - 1)

    # Calculate gamma_1 and gamma_2
    gamma_1 = (term1 - term2) / denominator
    gamma_2 = (term1 + term2) / denominator

    return gamma_1, gamma_2

if __name__ == '__main__':
    # Set the condition number kappa. You can change this value.
    kappa = 10.0

    # Get the optimal step sizes
    gamma_1, gamma_2 = find_optimal_steps(kappa)

    print(f"For a function with condition number kappa = {kappa}:")
    print(f"The best choice for the pair (gamma_1, gamma_2) is:")
    print(f"gamma_1 = {gamma_1}")
    print(f"gamma_2 = {gamma_2}")

    # The problem asks to output each number in the final equation.
    # Here are the components of the calculation for clarity.
    denominator = kappa**2 + 6 * kappa + 1
    term1_num = 4 * (kappa + 1)
    term2_num = 2 * math.sqrt(2) * (kappa - 1)

    print("\n--- Details of the final equations ---")
    print(f"The final equation for gamma_1 is: ({term1_num:.4f} - {term2_num:.4f}) / {denominator:.4f}")
    print(f"The final equation for gamma_2 is: ({term1_num:.4f} + {term2_num:.4f}) / {denominator:.4f}")
