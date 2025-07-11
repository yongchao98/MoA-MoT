def solve():
    """
    This function calculates and prints the formula for P(n).
    """
    
    # The coefficient of the n^-2 term in P(n)
    c2_num = 1
    c2_den = 720
    
    # The coefficient of the n^-4 term is derived from the expansion of exp(a2*n^-2 + a4*n^-4),
    # which is 1 + a2*n^-2 + (a4 + a2^2/2)*n^-4 + ...
    # a2 = 1/720
    # a4 = -1/5040
    # The coefficient is -1/5040 + (1/720)^2 / 2
    c4_num = -1433
    c4_den = 7257600
    
    # Print the formula for P(n)
    # The problem asks to output each number in the final equation.
    print(f"P(n) = 1 + {c2_num}/({c2_den} * n**2) - {abs(c4_num)}/({c4_den} * n**4)")

solve()