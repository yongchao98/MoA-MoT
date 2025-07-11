def solve_feynman_diagrams():
    """
    Calculates a(n), the number of non-vanishing Feynman diagrams of order 2n.
    The user asks for a(3).
    """
    n = 3
    
    # The formula is a(n) = (2n - 1)!!
    k = 2 * n - 1
    
    result = 1
    factors = []
    
    # Calculate the double factorial and store the factors
    for i in range(k, 0, -2):
        result *= i
        factors.append(str(i))
        
    # Create the string representation of the equation
    equation_str = " * ".join(factors)
    
    # Print the final result showing the full equation
    print(f"For n = {n}, the number of diagrams a({n}) is calculated as ({2*n - 1})!!.")
    print(f"a({n}) = {equation_str} = {result}")

solve_feynman_diagrams()