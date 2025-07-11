import math

def solve_expression():
    """
    Calculates the value of the simplified mathematical expression.
    """
    # The problem asks to calculate the value of:
    # ( (2/3)*l1*exp(l1/2) - (1/3)*exp(l1/2) ) * x2(0) - (2/3)*l2*exp(l2/2) - (10/3)*exp(l1/2)
    #
    # The value of x2(0) is not given. We assume the problem is constructed
    # such that the terms involving x2(0) cancel out with other terms.
    #
    # Let's group terms:
    # exp(l1/2) * [ ((2/3)*l1 - 1/3)*x2(0) - 10/3 ] - (2/3)*l2*exp(l2/2)
    #
    # Given l1 = (1 + sqrt(17))/2, we have (2*l1 - 1) = sqrt(17).
    # The term in the brackets becomes: [ (sqrt(17)/3)*x2(0) - 10/3 ]
    #
    # This term becomes zero if x2(0) = 10 / sqrt(17).
    # With this assumption, the expression simplifies to: - (2/3) * lambda_2 * exp(lambda_2 / 2)
    
    # Define constants
    sqrt17 = math.sqrt(17)
    coeff = -2/3

    # Define lambda_2
    lambda_2 = (1 - sqrt17) / 2

    # Calculate the exponential term
    exp_term = math.exp(lambda_2 / 2)

    # Calculate the final result
    result = coeff * lambda_2 * exp_term

    print("Based on a simplifying assumption, the final expression to calculate is: (-2/3) * lambda_2 * exp(lambda_2 / 2)")
    print("\nCalculating each number in the final equation:")
    print(f"The coefficient is -2/3 = {coeff:.6f}")
    print(f"lambda_2 = (1 - sqrt(17))/2 = {lambda_2:.6f}")
    print(f"The exponential term is exp(lambda_2 / 2) = {exp_term:.6f}")
    
    print("\nThe final equation is:")
    # The instruction is to "output each number in the final equation"
    print(f"({coeff}) * ({lambda_2}) * ({exp_term}) = {result}")

solve_expression()