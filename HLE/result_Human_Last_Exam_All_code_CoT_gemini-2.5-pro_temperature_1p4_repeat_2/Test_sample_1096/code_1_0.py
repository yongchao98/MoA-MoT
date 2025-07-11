import math

def calculate_optimal_steps(kappa):
    """
    Calculates the optimal two-step gradient descent learning rates.

    Args:
        kappa (float): The condition number (M/mu), assumed to be >= 1.

    Returns:
        tuple: A pair of floats (gamma1, gamma2) representing the optimal step sizes.
    """
    if kappa < 1:
        raise ValueError("Condition number kappa must be >= 1")
    if kappa == 1:
        # For kappa=1, the problem is trivial. The formula gives 1/1 = 1.
        # Let's see the limit: kappa -> 1. Denominator becomes 2*sqrt(2).
        # gamma1 = gamma2 = 1.
        return (1.0, 1.0)
        
    sqrt2 = math.sqrt(2)

    # These are the formulas derived from the Chebyshev polynomial analysis
    # for mu=1 and M=kappa.
    # The order of gamma1 and gamma2 can be swapped.
    
    # Formula for gamma1
    gamma1_num = 2 * sqrt2
    gamma1_den = kappa * (sqrt2 + 1) + (sqrt2 - 1)
    gamma1 = gamma1_num / gamma1_den
    
    # Formula for gamma2
    gamma2_num = 2 * sqrt2
    gamma2_den = kappa * (sqrt2 - 1) + (sqrt2 + 1)
    gamma2 = gamma2_num / gamma2_den

    return gamma1, gamma2

def main():
    """
    Main function to demonstrate the calculation for a sample kappa.
    """
    # Let's use a sample value for kappa, e.g., kappa = 10
    kappa = 10.0
    
    gamma1, gamma2 = calculate_optimal_steps(kappa)
    
    print(f"For a condition number kappa = {kappa}:")
    print("-" * 35)
    
    # The problem asks to output the numbers in the final equation.
    # Here are the equations for gamma_1 and gamma_2 with kappa substituted.
    print("The formulas for the optimal step sizes are:")
    sqrt2_str = "sqrt(2)"
    gamma1_formula = f"gamma_1 = (2 * {sqrt2_str}) / ({kappa} * ({sqrt2_str} + 1) + ({sqrt2_str} - 1))"
    gamma2_formula = f"gamma_2 = (2 * {sqrt2_str}) / ({kappa} * ({sqrt2_str} - 1) + ({sqrt2_str} + 1))"
    
    print(gamma1_formula)
    print(gamma2_formula)
    print("\nNumerically, these values are:")
    print(f"gamma_1 = {gamma1}")
    print(f"gamma_2 = {gamma2}")
    
    print("\nNote: The pair (gamma_1, gamma_2) is the set of optimal step sizes.")
    print("The order can be swapped between the two steps.")
    print("\nThe provided definition S = sqrt(M^2+(M-1)^2) was not used as it does not appear in the standard derivation for this problem and may be extraneous.")


if __name__ == "__main__":
    main()
