import math

def calculate_optimal_gammas(kappa):
    """
    Calculates the optimal two-step gradient descent step sizes (gamma_1, gamma_2)
    for a given condition number kappa.

    Args:
        kappa (float): The condition number of the function (M/mu).
                       Must be greater than 1.
    """
    if kappa <= 1:
        print("Error: The condition number kappa must be greater than 1.")
        return

    # Denominator from the derived formula
    denominator = kappa**2 + 6*kappa + 1

    # Numerator terms
    term1 = 4 * (kappa + 1)
    term2 = 2 * math.sqrt(2) * (kappa - 1)

    # Calculate gamma_1 and gamma_2
    gamma1_num = term1 + term2
    gamma2_num = term1 - term2

    gamma1 = gamma1_num / denominator
    gamma2 = gamma2_num / denominator

    # --- Outputting the results as requested ---
    print(f"For a function with condition number kappa = {kappa}:")
    print("-" * 40)
    
    print("The optimal step sizes (gamma_1, gamma_2) are determined by the formulas:")
    print("gamma_1 = (4*(kappa+1) + 2*sqrt(2)*(kappa-1)) / (kappa^2 + 6*kappa + 1)")
    print("gamma_2 = (4*(kappa+1) - 2*sqrt(2)*(kappa-1)) / (kappa^2 + 6*kappa + 1)")
    print("")

    # Print the detailed calculation for gamma_1
    print("Calculation for gamma_1:")
    print(f"gamma_1 = (4*({kappa}+1) + 2*{math.sqrt(2):.4f}*({kappa}-1)) / ({kappa}^2 + 6*{kappa} + 1)")
    print(f"gamma_1 = ({term1} + {term2:.4f}) / ({denominator})")
    print(f"gamma_1 = {gamma1_num:.4f} / {denominator}")
    print(f"gamma_1 = {gamma1:.6f}\n")
    
    # Print the detailed calculation for gamma_2
    print("Calculation for gamma_2:")
    print(f"gamma_2 = (4*({kappa}+1) - 2*{math.sqrt(2):.4f}*({kappa}-1)) / ({kappa}^2 + 6*{kappa} + 1)")
    print(f"gamma_2 = ({term1} - {term2:.4f}) / ({denominator})")
    print(f"gamma_2 = {gamma2_num:.4f} / {denominator}")
    print(f"gamma_2 = {gamma2:.6f}\n")

if __name__ == '__main__':
    # You can change this value to test with different condition numbers
    # As per the problem statement, we use M for the condition number kappa.
    M_condition_number = 10.0
    calculate_optimal_gammas(M_condition_number)
