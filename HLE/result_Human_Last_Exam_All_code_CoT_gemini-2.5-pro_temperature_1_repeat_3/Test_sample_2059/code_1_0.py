def solve():
    """
    Calculates the sum of the squares of the coefficients of the given polynomial expansion.
    """
    # The sum of squares is the constant term of P(x)P(1/x).
    # Let A_j be the constant term of the product for i from 0 to j.
    # We derived the recurrence relation A_j = 6*A_{j-1} - 2*A_{j-2}.
    # The base cases are A_{-1} = 1 (for an empty product) and A_0 = 4.
    
    # We store the computed values of the sequence A in a dictionary.
    a = {}
    a[-1] = 1
    a[0] = 4
    
    # We need to compute up to A_19.
    # The loop runs for j from 1 to 19.
    for j in range(1, 20):
        a[j] = 6 * a[j-1] - 2 * a[j-2]
        
    # The final answer is A_19.
    a19 = a[19]
    a18 = a[18]
    a17 = a[17]
    
    # The problem asks to output each number in the final equation.
    # We will show the final step of the recurrence relation as the equation.
    print(f"{a19} = 6 * {a18} - 2 * {a17}")

solve()