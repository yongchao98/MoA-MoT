import math

def solve_expression():
    """
    This function solves the problem as described.
    1. It determines x2(0) based on a steady-state analysis.
    2. It calculates the value of the given expression.
    3. It prints the full equation with numerical values and the final result.
    """
    
    # Values of lambda1 and lambda2 given in the problem
    lambda1 = (1 + math.sqrt(17)) / 2
    lambda2 = (1 - math.sqrt(17)) / 2
    
    # From the step-by-step analysis, we found that a constant equilibrium
    # solution x1(t)=6, x2(t)=0 satisfies all the system's conditions.
    # Therefore, we conclude that x2(0) = 0.
    x2_0 = 0
    
    # The expression to evaluate is:
    # (2/3 * lambda1 * exp(lambda1/2) - 1/3 * exp(lambda1/2)) * x2_0 
    # - 2/3 * lambda2 * exp(lambda2/2) - 10/3 * exp(lambda1/2)
    
    # Coefficients from the expression
    c1 = 2/3
    c2 = -1/3
    c3 = -2/3
    c4 = -10/3
    t = 1/2

    # Calculate intermediate values
    exp_l1_t = math.exp(lambda1 * t)
    exp_l2_t = math.exp(lambda2 * t)
    
    # Calculate each term of the expression
    term1 = (c1 * lambda1 * exp_l1_t + c2 * exp_l1_t) * x2_0
    term2 = c3 * lambda2 * exp_l2_t
    term3 = c4 * exp_l1_t
    
    result = term1 + term2 + term3
    
    # Print the equation with all the numbers, as requested
    print("The final equation is:")
    print(f"({c1:.4f} * {lambda1:.4f} * exp({lambda1:.4f}/2) {c2:.4f} * exp({lambda1:.4f}/2)) * {x2_0} "
          f"+ ({c3:.4f} * {lambda2:.4f} * exp({lambda2:.4f}/2)) "
          f"+ ({c4:.4f} * exp({lambda1:.4f}/2))")
          
    # With x2(0) = 0, the equation simplifies to:
    print("\nSince x2(0) = 0, the equation simplifies to:")
    print(f"({c3:.4f} * {lambda2:.4f} * exp({lambda2:.4f}/2)) + ({c4:.4f} * exp({lambda1:.4f}/2))")
    
    print(f"\nValue = {term2} + {term3}")
    print(f"Result = {result}")
    
    return result

if __name__ == '__main__':
    final_value = solve_expression()
    # The final output is requested in a specific format
    # print(f"<<<{final_value}>>>")