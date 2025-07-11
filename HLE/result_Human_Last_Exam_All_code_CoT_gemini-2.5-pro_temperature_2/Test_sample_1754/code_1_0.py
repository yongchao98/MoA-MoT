import math

def solve_expression():
    """
    Calculates the value of the given mathematical expression based on the provided context.
    """
    
    # Define the constants lambda_1 and lambda_2
    lambda1 = (1 + math.sqrt(17)) / 2
    lambda2 = (1 - math.sqrt(17)) / 2
    
    # The full expression is:
    # val = ( (2/3)*lambda1*e^(lambda1/2) - (1/3)*e^(lambda1/2) ) * x2(0) 
    #       - (2/3)*lambda2*e^(lambda2/2) 
    #       - (10/3)*e^(lambda1/2)
    
    # From the context of the controllability problem, we assume the system starts
    # at the origin, so x(0) = 0, which means x2(0) = 0.
    x2_0 = 0
    
    # With x2(0) = 0, the expression simplifies.
    # simplified_val = - (2/3)*lambda2*e^(lambda2/2) - (10/3)*e^(lambda1/2)
    
    # Calculate the two terms of the simplified expression
    term1_val = -(2/3) * lambda2 * math.exp(lambda2 / 2)
    term2_val = -(10/3) * math.exp(lambda1 / 2)
    
    # The final value is the sum of these two terms
    final_value = term1_val + term2_val
    
    # As requested, output each number in the final equation.
    # The equation is: final_value = term1 + term2
    print("Assuming x₂(0) = 0, the expression simplifies to:")
    print("-(2/3) * λ₂ * e^(λ₂/2) - (10/3) * e^(λ₁/2)")
    print("\nThe equation with the calculated values for each term is:")
    print(f"{term1_val:.7f} + ({term2_val:.7f}) = {final_value:.7f}")

solve_expression()