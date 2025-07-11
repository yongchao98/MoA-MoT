import math

def solve_problem():
    """
    Calculates the final value based on the derived formula.
    
    The formula is: (n * (n + 1)) / (2 * pi * (n - 1)) + 10^15
    """
    n = 4048.0
    val_10_15 = 10.0**15

    # Components of the first term in the expression
    numerator = n * (n + 1.0)
    denominator = 2.0 * math.pi * (n - 1.0)
    
    # Calculate the final result
    result = numerator / denominator + val_10_15
    
    # Output the components of the calculation as requested
    print("The final expression to calculate is: (n * (n + 1)) / (2 * pi * (n - 1)) + 10^15")
    print("\nWhere the values are:")
    print(f"n = {int(n)}")
    print(f"n + 1 = {int(n + 1.0)}")
    print(f"2 = 2")
    print(f"pi = {math.pi}")
    print(f"n - 1 = {int(n - 1.0)}")
    print(f"10^15 = {int(val_10_15)}")
    
    # Print the equation with values plugged in
    print(f"\nCalculation:")
    print(f"({int(n)} * {int(n + 1.0)}) / (2 * {math.pi} * {int(n - 1.0)}) + {int(val_10_15)}")
    
    # Print the final numerical answer
    print(f"\nResult: {result}")
    
    return result

# Execute the function to get the final answer
final_answer = solve_problem()