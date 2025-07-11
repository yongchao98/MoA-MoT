def solve_feynman_diagrams():
    """
    Calculates a(3), the number of non-vanishing Feynman diagrams of order 6
    in quantum electrodynamics.
    The formula is a(n) = (2n-1)!!.
    """
    n = 3
    
    # The number to calculate the double factorial of
    k = 2 * n - 1
    
    result = 1
    factors = []
    
    # Calculate the double factorial (k!!)
    for i in range(k, 0, -2):
        result *= i
        factors.append(str(i))
        
    # Format the equation string for printing
    equation_str = " * ".join(factors)
    
    print(f"The number of non-vanishing Feynman diagrams, a(n), of order 2n is given by the formula a(n) = (2n-1)!!.")
    print(f"For n = {n}, the order is 2n = {2*n}, and we need to calculate a({n}).")
    print(f"a({n}) = (2*{n} - 1)!! = {k}!!")
    print(f"The final equation is: {equation_str} = {result}")

solve_feynman_diagrams()