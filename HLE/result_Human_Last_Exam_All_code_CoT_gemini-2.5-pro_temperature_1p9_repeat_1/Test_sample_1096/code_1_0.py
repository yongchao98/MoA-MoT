import math

def calculate_optimal_stepsizes(kappa):
    """
    Calculates the optimal two-step gradient descent learning rates (gamma1, gamma2).

    Args:
        kappa (float): The condition number M/mu of the function. Assumes mu=1, M=kappa.
                       Must be greater than 1.
    """
    if not isinstance(kappa, (int, float)) or kappa <= 1:
        print("Error: kappa must be a number greater than 1.")
        return

    # Denominator from the derived formula
    # D = kappa^2 + 6*kappa + 1
    denominator = kappa**2 + 6 * kappa + 1

    # Calculate the common part of the numerators
    term1 = 4 * (kappa + 1)
    term2 = 2 * math.sqrt(2) * (kappa - 1)

    # Calculate gamma_1 and gamma_2 using the derived formulas
    # The order of gamma_1 and gamma_2 can be swapped without affecting the final result.
    gamma1 = (term1 - term2) / denominator
    gamma2 = (term1 + term2) / denominator

    print(f"For a condition number kappa = {kappa}:")
    print(f"The optimal learning rate pair (gamma_1, gamma_2) is:")
    print(f"gamma_1 = {gamma1}")
    print(f"gamma_2 = {gamma2}")


if __name__ == '__main__':
    # You can change this value to test different condition numbers.
    # The condition number M is denoted as kappa, assuming mu=1.
    kappa = 10.0
    calculate_optimal_stepsizes(kappa)
    # Example for another condition number
    # kappa_2 = 100.0
    # calculate_optimal_stepsizes(kappa_2)