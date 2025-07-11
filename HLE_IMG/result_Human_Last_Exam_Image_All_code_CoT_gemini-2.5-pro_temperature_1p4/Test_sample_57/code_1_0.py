def print_tiling_formula():
    """
    This function prints the derived formula for the number of ways
    to tile the L-shaped figure. The formula is expressed in terms of
    Fibonacci numbers.
    """
    # The formula is F_n^2 + F_{n-1}^2
    # The variable is n. The numbers in the formula are 2 and 1.
    
    variable = "n"
    exponent = 2
    operator = "+"
    one = 1
    
    print("The number of ways to tile the shape for a given n is:")
    
    # We print the formula piece by piece to highlight the numbers involved.
    # Final form: F_n^2 + F_{n-1}^2
    
    print(f"F_{variable}^{exponent} {operator} F_({variable}-{one})^{exponent}")

print_tiling_formula()