def solve():
    """
    This function provides the formula for the correction factor P(n).
    """
    
    # Coefficients for the terms in P(n)
    c0 = 1
    
    # Coefficient for the 1/n^2 term
    c1_num = 1
    c1_den = 720
    
    # Coefficient for the 1/n^4 term
    c2_num = -1433
    c2_den = 7257600

    # Construct and print the formula for P(n)
    # The requirement is to output each number in the final equation.
    print(f"P(n) = {c0} + {c1_num}/({c1_den}*n^2) - {abs(c2_num)}/({c2_den}*n^4)")

solve()