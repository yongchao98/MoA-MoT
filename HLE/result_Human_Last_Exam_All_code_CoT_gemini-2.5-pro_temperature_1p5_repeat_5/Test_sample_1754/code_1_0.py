import math

def solve_expression():
    """
    Calculates the value of the given mathematical expression.
    """
    # Given constants
    lambda_1 = (1 + math.sqrt(17)) / 2
    lambda_2 = (1 - math.sqrt(17)) / 2

    # In controllability problems, it's a common task to analyze reachability from the origin.
    # Therefore, we assume the initial condition x(0) = 0, which means x_2(0) = 0.
    x2_0 = 0

    # Pre-calculate the exponential terms for clarity
    exp_lambda1_half = math.exp(lambda_1 / 2)
    exp_lambda2_half = math.exp(lambda_2 / 2)
    
    # The full expression is:
    # (2/3 * lambda_1 * exp(lambda_1/2) - 1/3 * exp(lambda_1/2)) * x_2(0)
    # - 2/3 * lambda_2 * exp(lambda_2/2) - 10/3 * exp(lambda_1/2)
    
    # Calculate each part of the expression
    term1_coeff = (2/3) * lambda_1 * exp_lambda1_half - (1/3) * exp_lambda1_half
    term1 = term1_coeff * x2_0
    
    term2 = -(2/3) * lambda_2 * exp_lambda2_half
    
    term3 = -(10/3) * exp_lambda1_half
    
    # The final result is the sum of these terms
    result = term1 + term2 + term3

    # Output the explanation and calculation steps
    print("Assuming the initial state x(0) = 0, which implies x_2(0) = 0.")
    print("The expression to calculate is:")
    print("(2/3 * lambda_1 * exp(lambda_1/2) - 1/3 * exp(lambda_1/2)) * x_2(0) - 2/3 * lambda_2 * exp(lambda_2/2) - 10/3 * exp(lambda_1/2)")
    
    print("\nSubstituting the values:")
    print(f"lambda_1 = {lambda_1:.5f}")
    print(f"lambda_2 = {lambda_2:.5f}")
    print(f"exp(lambda_1/2) = {exp_lambda1_half:.5f}")
    print(f"exp(lambda_2/2) = {exp_lambda2_half:.5f}")
    print(f"x_2(0) = {x2_0}")

    print("\nThe equation with numbers becomes:")
    # Print the equation with all numbers substituted, as requested.
    print(f"({2/3:.5f} * {lambda_1:.5f} * {exp_lambda1_half:.5f} - {1/3:.5f} * {exp_lambda1_half:.5f}) * {x2_0} "
          f"- {2/3:.5f} * {lambda_2:.5f} * {exp_lambda2_half:.5f} "
          f"- {10/3:.5f} * {exp_lambda1_half:.5f}")

    # Print the calculated values for each major term
    print(f"\n= ({term1_coeff:.5f}) * {x2_0} + ({term2:.5f}) + ({term3:.5f})")
    
    # Print the final result
    print(f"\n= {result:.5f}")

solve_expression()