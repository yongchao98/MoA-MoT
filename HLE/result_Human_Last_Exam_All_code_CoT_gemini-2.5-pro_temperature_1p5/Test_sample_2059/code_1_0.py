def solve():
    """
    Solves the problem by calculating the sum of the squares of the coefficients
    of the given polynomial product.
    """
    n = 19
    
    # We need to compute A_n and B_n such that (3 + sqrt(7))^n = A_n + B_n * sqrt(7).
    # We can find a recurrence relation for A_n and B_n.
    # A_n = 3 * A_{n-1} + 7 * B_{n-1}
    # B_n = A_{n-1} + 3 * B_{n-1}
    # The base case is n=0: (3 + sqrt(7))^0 = 1 = 1 + 0*sqrt(7), so A_0=1, B_0=0.
    
    a = 1  # A_0
    b = 0  # B_0
    
    for _ in range(n):
        a_new = 3 * a + 7 * b
        b_new = a + 3 * b
        a = a_new
        b = b_new
        
    # The value we are looking for is V_19, which is derived from the recurrence relation.
    # The final formula for V_n is 4*A_n + 10*B_n.
    # So we need to calculate 4*A_19 + 10*B_19.
    
    a_19 = a
    b_19 = b
    
    result = 4 * a_19 + 10 * b_19
    
    print(f"The calculation is based on the formula: 4 * A_{n} + 10 * B_{n}")
    print(f"where A_{n} + B_{n} * sqrt(7) = (3 + sqrt(7))^{n}, and n = {n}.")
    print(f"A_{n} = {a_19}")
    print(f"B_{n} = {b_19}")
    print(f"Result = 4 * {a_19} + 10 * {b_19} = {result}")

solve()