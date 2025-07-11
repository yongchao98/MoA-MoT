def solve():
    """
    Calculates the complexity representation (a,b,c) for the two regimes.

    The query complexity for a given sorting problem is represented in the form
    Theta(sqrt(N^a * (log N)^b * (log log N)^c)).

    If the derived complexity is Theta(N^p * (log N)^q * (log log N)^r),
    then by squaring it and comparing with the expression inside the square root, we get:
    a = 2*p
    b = 2*q
    c = 2*r

    Our analysis shows that for both specified regimes, the optimal sorting strategy
    results in a query complexity of Theta(N * log N).
    """

    # For Theta(N * log N) = Theta(N^1 * (log N)^1 * (log log N)^0),
    # the powers are p=1, q=1, r=0.
    p = 1
    q = 1
    r = 0

    # Calculate a, b, c for the first regime
    a1 = 2 * p
    b1 = 2 * q
    c1 = 2 * r
    result1 = (a1, b1, c1)

    # Calculate a, b, c for the second regime
    # The complexity is the same
    a2 = 2 * p
    b2 = 2 * q
    c2 = 2 * r
    result2 = (a2, b2, c2)
    
    # Print the final result in the requested format "(a,b,c),(d,e,f)"
    # We output each number of the final tuples
    print(f"({result1[0]},{result1[1]},{result1[2]}),({result2[0]},{result2[1]},{result2[2]})")

solve()