import math

def solve_imitation_learning_bound():
    """
    This function calculates and prints the tightest upper bound for the
    performance difference J(pi^*) - J(pi_hat) based on the problem statement.
    """

    # The final expression for the upper bound is H^2 * |A| * (1 - exp(-lambda)).
    # We will print this formula as a string.
    
    bound_formula = "H**2 * |A| * (1 - exp(-lambda))"

    print("The expression for the tightest upper bound of J(pi^*) - J(pi_hat) is:")
    print(bound_formula)
    
    # The prompt requests to output each number in the final equation.
    # The numbers in the derived formula are the exponent 2 for H, the constant 1,
    # and the constant -1 multiplying lambda in the exponent.
    
    print("\nThe numbers present in this final equation are:")
    print(f"The exponent on H: 2")
    print(f"The constant in the term (1 - ...): 1")
    print(f"The constant factor of lambda in the exponent: -1")

solve_imitation_learning_bound()