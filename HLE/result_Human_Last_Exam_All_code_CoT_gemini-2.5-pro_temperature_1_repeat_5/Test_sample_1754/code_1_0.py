import math

def solve_expression():
    """
    Calculates the value of the given mathematical expression.
    """
    # Given lambda values
    lambda1 = (1 + math.sqrt(17)) / 2
    lambda2 = (1 - math.sqrt(17)) / 2
    
    # In controllability problems, it's common to assume the system starts
    # from the origin, hence x(0) = 0. This implies x2(0) = 0.
    x2_0 = 0

    print("The expression to calculate is:")
    print("(2/3 * lambda1 * exp(lambda1/2) - 1/3 * exp(lambda1/2)) * x2(0) - 2/3 * lambda2 * exp(lambda2/2) - 10/3 * exp(lambda1/2)\n")
    
    print(f"Assuming the system starts from the origin, we have x2(0) = {x2_0}.")
    print("The expression simplifies to: - (2/3) * lambda2 * exp(lambda2/2) - (10/3) * exp(lambda1/2)\n")

    # Calculate the components of the expression
    term1_val = (-2/3) * lambda2 * math.exp(lambda2 / 2)
    term2_val = (-10/3) * math.exp(lambda1 / 2)
    
    total_value = term1_val + term2_val

    print("--- Calculation Steps ---")
    print(f"lambda1 = {lambda1}")
    print(f"lambda2 = {lambda2}")
    print(f"exp(lambda1/2) = {math.exp(lambda1 / 2)}")
    print(f"exp(lambda2/2) = {math.exp(lambda2 / 2)}")
    
    print("\n--- Final Equation ---")
    # Print each number in the final equation
    print(f"Expression: ({(-2/3):.6f}) * ({lambda2:.6f}) * ({math.exp(lambda2 / 2):.6f}) + ({(-10/3):.6f}) * ({math.exp(lambda1 / 2):.6f})")
    print(f"Calculation: {term1_val:.6f} + {term2_val:.6f}")
    print(f"Result: {total_value:.6f}")
    
    print(f"\nThe final value is: {total_value}")
    
    return total_value

# Execute the function
final_answer = solve_expression()