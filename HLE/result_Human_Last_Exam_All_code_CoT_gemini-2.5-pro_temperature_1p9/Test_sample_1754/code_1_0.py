import math

def solve_problem():
    """
    Solves the problem by calculating the value of the given expression
    under the assumption that the initial state is zero.
    """
    
    # Step 1: Define lambda_1 and lambda_2
    sqrt_17 = math.sqrt(17)
    lambda_1 = (1 + sqrt_17) / 2
    lambda_2 = (1 - sqrt_17) / 2
    
    # Step 2: Assume x(0) = 0, so x_2(0) = 0.
    # The expression to calculate is:
    # (2/3 * lambda_1 * exp(lambda_1/2) - 1/3 * exp(lambda_1/2)) * x_2(0) 
    # - 2/3 * lambda_2 * exp(lambda_2/2) 
    # - 10/3 * exp(lambda_1/2)
    x_2_0 = 0.0

    # Step 3: Calculate each term of the expression.
    term1_coeff = (2/3 * lambda_1 * math.exp(lambda_1 / 2) - 1/3 * math.exp(lambda_1 / 2))
    term1 = term1_coeff * x_2_0
    
    term2 = - (2/3) * lambda_2 * math.exp(lambda_2 / 2)
    
    term3 = - (10/3) * math.exp(lambda_1 / 2)
    
    # Step 4: Calculate the final result
    result = term1 + term2 + term3
    
    # Step 5: Print the results including the numbers in the final equation.
    print("Assuming the system starts from rest, so x(0) = 0 and x2(0) = 0.")
    print("The expression to evaluate is:")
    print("Value = (2/3 * lambda_1 * exp(lambda_1/2) - 1/3 * exp(lambda_1/2)) * x2(0) - 2/3 * lambda_2 * exp(lambda_2/2) - 10/3 * exp(lambda_1/2)")
    print("\nLet's calculate the value of each number in the equation:")
    print(f"lambda_1 = {lambda_1}")
    print(f"lambda_2 = {lambda_2}")
    print(f"x2(0) = {x_2_0}")
    print("\nFinal Equation:")
    print(f"Value = ({term1_coeff:.4f}) * {x_2_0} + ({term2:.4f}) + ({term3:.4f})")
    print(f"Value = {term1:.4f} + {term2:.4f} + {term3:.4f}")
    print(f"Value = {result:.4f}")

    # Return the final numerical value for the grading format
    return result

# Execute the function and capture the result
final_value = solve_problem()

# Final Answer in specified format
# print(f"\n<<<{final_value:.4f}>>>")