import sympy

def solve_l_n():
    """
    This function calculates the symbolic expression for l(n) using sympy.
    """
    # Define n as a symbol, for n >= 5
    n = sympy.Symbol('n', integer=True, positive=True)

    # Step 1: Define matrices M and A (which is P^{-1})
    # M has diagonal a and off-diagonal b
    a = sympy.sqrt(1 - (n - 1) / n**2)
    b = 1/n

    # A is tridiagonal with 2 on the diagonal and 1 on the super/sub-diagonals.
    
    # Step 2: Calculate the projection coefficients d_j
    # d_j = (MA)_{jj}
    # d_1 = 2a+b
    # d_n = 2a+b
    # d_i = 2a+2b for i in {2, ..., n-1}
    d1 = 2*a + b
    dn = 2*a + b
    di = 2*a + 2*b

    # Step 3: Calculate the sum S1 for the first row adjustment
    # S1 = sum_j M_1j * d_j
    # S1 = M_11*d1 + M_12*d2 + ... + M_1n*dn
    # M_11=a, M_1j=b for j>1
    # The sum is a*d1 + b*d2 + ... + b*d(n-1) + b*dn
    # which is a*d1 + b * (n-2)*di + b*dn
    S1 = a*d1 + (n-2)*b*di + b*dn
    
    # By symmetry, S_n = S_1

    # Step 4: Calculate l(n)
    # Sum of A's first row is 3. Sum of A's last row is 3.
    # l(n) = (sum of row 1 of B) + (sum of row n of B)
    #      = (sum of row 1 of A - S1) + (sum of row n of A - Sn)
    #      = (3 - S1) + (3 - Sn) = 6 - 2*S1
    l_n = 6 - 2*S1

    # Simplify the expression
    l_n_simplified = sympy.simplify(l_n)

    # Print the result in a readable format, showing the numbers in the equation
    # The simplified expression is (2*(n**2 + 1) - 2*(2*n - 1)*sqrt(n**2 - n + 1))/n**2
    num_part1 = 2 * (n**2 + 1)
    num_part2_coeff = 2 * (2*n - 1)
    num_part2_sqrt = sympy.sqrt(n**2 - n + 1)
    denom = n**2
    
    print("The exact value of l(n) is given by the expression:")
    print(f"l(n) = ( {num_part1} - ({num_part2_coeff})*sqrt({n**2 - n + 1}) ) / ( {denom} )")

    # To satisfy the prompt's instruction: "output each number in the final equation"
    # We will print the coefficients and constants from the final expanded form:
    # l(n) = (2*n**2 + 2 - (4*n - 2)*sqrt(n**2 - n + 1)) / n**2
    print("\nNumbers in the numerator expression '2*n**2 + 2 - (4*n - 2)*sqrt(n**2 - n + 1)':")
    print(f"Coefficients and constants are: 2, 2, 2, 4, 2, 2, 1, 1")
    # Explanation: 2,2 (from 2n^2+2), 4,2 (from 4n-2), 2,1,1 (from n^2-n+1)

solve_l_n()