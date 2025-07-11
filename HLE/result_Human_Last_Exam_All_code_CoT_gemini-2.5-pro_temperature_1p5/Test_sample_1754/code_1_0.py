import math

def solve_problem():
    """
    Calculates the value of the given expression based on the problem statement.

    The problem is a controllability problem where the initial state x(0) is not
    explicitly given. For such problems, a standard approach is to assume the system
    starts from the origin, i.e., x(0) = 0. This makes x_2(0) = 0.

    With x_2(0) = 0, the expression to evaluate simplifies to:
    - (2/3) * lambda_2 * exp(lambda_2 / 2) - (10/3) * exp(lambda_1 / 2)
    """

    # Define the eigenvalues
    lambda_1 = (1 + math.sqrt(17)) / 2
    lambda_2 = (1 - math.sqrt(17)) / 2
    
    # Assume x_2(0) = 0
    x2_0 = 0

    # Calculate each term of the original expression
    term1_coeff_val = (2/3 * lambda_1 * math.exp(lambda_1 / 2) - 1/3 * math.exp(lambda_1 / 2))
    term1 = term1_coeff_val * x2_0
    
    term2_val = -2/3 * lambda_2 * math.exp(lambda_2 / 2)
    
    term3_val = -10/3 * math.exp(lambda_1 / 2)

    # Calculate the final result
    result = term1 + term2_val + term3_val

    # Print the equation with the calculated values
    print(f"Assuming x_2(0) = 0, the expression is:")
    print(f"({term1_coeff_val:.4f}) * ({x2_0}) + ({term2_val:.4f}) + ({term3_val:.4f})")
    
    # Print the step-by-step calculation
    print(f"= {term1:.4f} + ({term2_val:.4f}) + ({term3_val:.4f})")
    
    final_eq = f"= {term2_val:.4f} - {abs(term3_val):.4f}"
    
    # Print the final result
    print(f"{final_eq} = {result:.4f}")
    
    # Returning the final numeric result for the answer tag
    return result

# Execute the function and capture the final numeric result
final_answer = solve_problem()