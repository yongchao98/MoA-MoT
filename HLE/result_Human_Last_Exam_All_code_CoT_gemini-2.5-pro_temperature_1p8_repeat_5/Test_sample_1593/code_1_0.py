def solve():
    """
    This function prints the calculated complexity values in abc-notation for the two regimes.
    The derivation is as follows:

    Regime 1: N = 2^sqrt(L), which implies L = (log N)^2.
    A natural algorithm is a Radix Sort with block size b = log N.
    - Number of stages = L/b = log N.
    - H-queries cost = N * (L/b) = N * log N.
    - C-queries cost = (L/b) * O(N log N) = O(N * (log N)^2).
    - Total Query Complexity Q = Theta(N * (log N)^2).
    - Q^2 = Theta(N^2 * (log N)^4).
    - This corresponds to (a,b,c) = (2,4,0).

    Regime 2: N = 2^((log L)^2), which implies L = 2^sqrt(log N).
    Here L is very large and not a simple polynomial in log N. A natural algorithm is to
    treat the strings as atomic (Full-String Hashing Sort).
    - H-queries cost = N (to hash each full string).
    - C-queries cost = O(N log N) (to sort the unique strings).
    - Total Query Complexity Q = Theta(N * log N).
    - Q^2 = Theta(N^2 * (log N)^2).
    - This corresponds to (a,b,c) = (2,2,0).
    """
    
    # Complexity for the first regime: (a, b, c) = (2, 4, 0)
    a1, b1, c1 = 2, 4, 0

    # Complexity for the second regime: (a, b, c) = (2, 2, 0)
    a2, b2, c2 = 2, 2, 0

    # Print the result in the specified format "(a,b,c),(a,b,c)"
    print(f"({a1},{b1},{c1}),({a2},{b2},{c2})")

solve()